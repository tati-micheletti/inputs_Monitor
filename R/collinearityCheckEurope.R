#' Resolve final predictors and assemble the model-ready table (Europe scale)
#'
#' For every species, always saves a correlation-structure corrplot over
#' all candidate bioclim variables (diagnostic, independent of the toggle
#' below). `speciesConfig_predictors.csv`'s table (`speciesPredictorTable`)
#' is the ONLY source of a species' candidate predictors -- no mode
#' selector. Optionally, `dropCollinearPredictors` prunes that species' own
#' listed predictors for collinearity via real block-CV selection
#' (`select07Blockcv()`); when FALSE (default), the table's list is used
#' exactly as given.
#'
#' Regardless of the toggle, the output table always has the same
#' guaranteed base columns (cell50x50, birdlife_scientific_name,
#' occurrence, x, y, foldID) plus whatever predictor columns were resolved.
#'
#' @param pooledData Named list (by species) of occurrence+bioclim data.frames.
#' @param blocksData Named list (by species) of blocks objects (real or
#'   mimicked) -- must have `$folds_ids` aligned to `pooledData[[sp]]` rows.
#' @param speciesPredictorTable Named list (species -> character vector).
#'   Sourced from `speciesConfig_predictors.csv` via
#'   `loadSpeciesPredictorConfig()`. A species missing here is a hard error --
#'   there's no fallback source for a species' predictors.
#' @param dropCollinearPredictors Logical, default FALSE. If TRUE, runs real
#'   block-CV collinearity selection (`select07Blockcv()`) over each
#'   species' own table-listed predictors, capped at 1 predictor per 10
#'   occurrences -- pruning what the table says to consider, never
#'   substituting a different candidate set. If FALSE, the table's list is
#'   used exactly as given.
#' @param corrplotDir Character. Directory to save correlation plots in.
#' @param threshold Numeric. Absolute correlation threshold, default 0.7.
#' @param univar Character. Initial univariate model form, default "gam".
#' @param spatialTermSpecies Named list (species -> TRUE/FALSE), or NULL.
#'   Species with `TRUE` get projected `x`/`y` coordinates added as an extra
#'   predictor on top of the table's list -- a spatial trend-surface term,
#'   NOT a formal random effect (see DECISIONS.md's 2026-09-26 entries). A
#'   species absent from this list, or set `FALSE`, never gets it.
#'   Deliberately opt-in per species -- can just as easily hurt a model
#'   (overfitting to historical geography, reduced transportability to
#'   future predictions, diluted variable-importance interpretation) as
#'   help it.
#' @param cachePath Character, or NULL (default). Directory for
#'   `reproducible::Cache()`'s per-species cache (see
#'   `resolveEuropeSpeciesPredictors()`) -- e.g. `cachePath(sim)`, a
#'   stable location shared across runs. NULL falls back to a temp
#'   directory, for standalone/test calls.
#' @return Named list (by species) with `data` (the final table) and
#'   `predictors` (character vector of predictor columns used -- includes
#'   `x`/`y` only for species opted into `spatialTermSpecies`).
collinearityCheckEurope <- function(pooledData, blocksData, speciesPredictorTable,
                                     dropCollinearPredictors = FALSE,
                                     spatialTermSpecies = NULL, corrplotDir,
                                     threshold = 0.7, univar = "gam",
                                     cachePath = NULL) {

  if (is.null(cachePath)) cachePath <- file.path(tempdir(), "birdMonitor_cache")
  dir.create(corrplotDir, recursive = TRUE, showWarnings = FALSE)
  bioVars <- bioclimPredictorColumns()
  result <- list()

  for (sp in names(pooledData)) {
    message("\n  Processing: ", sp)
    spPa <- pooledData[[sp]]
    spClean <- gsub(" ", "_", sp)

    missing <- setdiff(bioVars, names(spPa))
    if (length(missing) > 0) {
      stop("Missing bioclim variables for ", sp, ": ", paste(missing, collapse = ", "))
    }

    corMat <- cor(spPa[, bioVars], method = "spearman")
    corrplotFile <- file.path(corrplotDir, paste0(spClean, "_EU.png"))
    grDevices::png(corrplotFile, width = 1600, height = 1600, res = 240)
    corrplot::corrplot.mixed(corMat, tl.pos = "lt", tl.cex = 0.6, number.cex = 0.5,
                              addCoefasPercent = TRUE)
    grDevices::dev.off()

    if (is.null(speciesPredictorTable[[sp]])) {
      stop(sp, ": no entry in speciesPredictorTable (speciesConfig_predictors.csv) -- ",
           "every species must be listed there, there is no fallback.")
    }

    spResult <- reproducible::Cache(
      resolveEuropeSpeciesPredictors, sp = sp, spPa = spPa, blocksSp = blocksData[[sp]],
      bioVars = bioVars, requestedPredictors = speciesPredictorTable[[sp]],
      dropCollinearPredictors = dropCollinearPredictors,
      spatialTerm = isTRUE(spatialTermSpecies[[sp]]), threshold = threshold, univar = univar,
      cachePath = cachePath, userTags = c("collinearityCheckEurope", spClean))

    result[[sp]] <- spResult
  }

  result
}

#' Resolve one species' final predictor set and assemble its model-ready table
#'
#' See `resolveHabitatSpeciesPredictors()` in `collinearityCheckGerHabitat.R`
#' for the full rationale (same pattern, Europe/climate scale).
#'
#' @param bioVars Character vector. Candidate bioclim predictor names.
#' @inheritParams resolveHabitatSpeciesPredictors
#' @return List with `data` (the final table) and `predictors`.
resolveEuropeSpeciesPredictors <- function(sp, spPa, blocksSp, bioVars, requestedPredictors,
                                            dropCollinearPredictors, spatialTerm,
                                            threshold, univar) {
  requested <- requestedPredictors
  requested <- requested[!is.na(requested) & nzchar(trimws(requested))]

  if (length(requested) == 0) {
    stop(sprintf(
      "Species '%s' has no climate predictors listed in speciesConfig_predictors.csv!",
      sp
    ))
  }

  missingCols <- setdiff(requested, bioVars)
  if (length(missingCols) > 0) {
    stop(sprintf(
      "\n[PREDICTOR ERROR] Species '%s' [CLIMATE scale]:\nRequested climate predictor(s) %s were NOT found in Europe dataset!\nAvailable climate predictors are: %s\n",
      sp,
      paste(dQuote(missingCols), collapse = ", "),
      paste(bioVars, collapse = ", ")
    ))
  }

  if (isTRUE(dropCollinearPredictors)) {
    nPres <- sum(spPa$occurrence == 1)
    nAbs <- sum(spPa$occurrence == 0)
    varSel <- select07Blockcv(X = spPa[, requested], y = spPa$occurrence, threshold = threshold,
                              univar = univar, spBlock = blocksSp, weights = rep(1, nrow(spPa)))
    occNum <- max(floor(min(nPres, nAbs) / 10), 1)
    predSel <- stats::na.omit(varSel$pred_sel[1:min(occNum, length(varSel$pred_sel))])
  } else {
    predSel <- requested
  }

  # Opt-in per species (see spatialTermSpecies docstring above) -- add
  # projected x/y coordinates as predictors, on top of the table's list.
  # Not a formal random effect (BRT/dismo::gbm.step() has no mixed-model
  # machinery) -- functionally, letting the tree split on location itself.
  if (isTRUE(spatialTerm)) {
    predSel <- c(as.character(predSel), "x", "y")
  }

  message("Predictors used (", length(predSel), "): ", paste(predSel, collapse = ", "))

  spPaOut <- spPa
  spPaOut$foldID <- blocksSp$folds_ids

  keepCols <- c("cell50x50", "birdlife_scientific_name", "occurrence", "x", "y", "foldID", predSel)
  list(data = spPaOut[, intersect(keepCols, names(spPaOut))],
       predictors = as.character(predSel))
}
