#' Resolve final predictors and assemble the model-ready table (German habitat scale)
#'
#' `speciesConfig_predictors.csv`'s table (`speciesPredictorTable`) is the ONLY
#' source of a species' candidate predictors -- no mode selector. Optionally,
#' `dropCollinearPredictors` prunes that species' own listed predictors for
#' collinearity via real block-CV selection (`select07Blockcv()`); when
#' FALSE (default), the table's list is used exactly as given. This
#' guarantees the same base output columns (AREA_NATIONAL_CODE, year,
#' latin_name, occurrence, x, y, foldID) regardless of the toggle.
#'
#' @param pooledData Named list (by species) of pooled occurrence+covariate
#'   data.frames (see `poolOccurrenceGerHabitat()`).
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
#'   `resolveHabitatSpeciesPredictors()`) -- e.g. `cachePath(sim)`, a
#'   stable location shared across runs. NULL falls back to a temp
#'   directory, for standalone/test calls.
#' @return Named list (by species) with `data` (the final table) and
#'   `predictors` (character vector of predictor columns used -- includes
#'   `x`/`y` only for species opted into `spatialTermSpecies`).
collinearityCheckGerHabitat <- function(pooledData, blocksData, speciesPredictorTable,
                                        dropCollinearPredictors = FALSE,
                                        spatialTermSpecies = NULL, corrplotDir,
                                        threshold = 0.7, univar = "gam",
                                        cachePath = NULL) {

  if (is.null(cachePath)) cachePath <- file.path(tempdir(), "birdMonitor_cache")
  dir.create(corrplotDir, recursive = TRUE, showWarnings = FALSE)
  allPredictors <- covariatePredictorColumns()
  result <- list()

  for (sp in names(pooledData)) {
    message("\n  -- ", sp, " --------------------------")
    spPa <- pooledData[[sp]]
    spClean <- gsub(" ", "_", sp)

    if (is.null(blocksData[[sp]])) {
      message("No blocks available for ", sp, " -- skipping")
      next
    }

    nPres <- sum(spPa$occurrence == 1)
    nAbs <- sum(spPa$occurrence == 0)
    message("Pooled records: ", nrow(spPa), " (", nPres, " pres / ", nAbs, " abs)")

    # Candidate environmental predictors present in data -- only computed
    # here for the diagnostic corrplot (a general diagnostic over every
    # candidate column, independent of this species' own predictor table);
    # resolveHabitatSpeciesPredictors() works from the table directly, not
    # this candidate list, so it stays self-contained under Cache().
    predCols <- intersect(allPredictors, names(spPa))
    allNaCols <- predCols[sapply(predCols, function(col) all(is.na(spPa[[col]])))]
    if (length(allNaCols) > 0) {
      message("Removing all-NA predictors: ", paste(allNaCols, collapse = ", "))
      predCols <- setdiff(predCols, allNaCols)
    }

    X <- spPa[, predCols, drop = FALSE]
    X <- X[, sapply(X, function(col) length(unique(col[!is.na(col)])) > 1), drop = FALSE]

    if (ncol(X) >= 2) {
      corrplotFile <- file.path(corrplotDir, paste0(spClean, "_habitat_pooled.png"))
      if (!file.exists(corrplotFile)) {
        grDevices::png(corrplotFile, width = 1200, height = 1200, res = 150)
        corrplot::corrplot.mixed(cor(X, method = "spearman", use = "complete.obs"),
                                 tl.pos = "lt", tl.cex = 0.6, number.cex = 0.4, addCoefasPercent = TRUE)
        grDevices::dev.off()
      }
    }

    if (is.null(speciesPredictorTable[[sp]])) {
      stop(sp, ": no entry in speciesPredictorTable (speciesConfig_predictors.csv) -- ",
           "every species must be listed there, there is no fallback.")
    }

    spResult <- reproducible::Cache(
      resolveHabitatSpeciesPredictors, sp = sp, spPa = spPa, blocksSp = blocksData[[sp]],
      requestedPredictors = speciesPredictorTable[[sp]],
      dropCollinearPredictors = dropCollinearPredictors,
      spatialTerm = isTRUE(spatialTermSpecies[[sp]]), threshold = threshold, univar = univar,
      cachePath = cachePath, userTags = c("collinearityCheckGerHabitat", spClean))

    result[[sp]] <- spResult
  }

  result
}

#' Resolve one species' final predictor set and assemble its model-ready table
#'
#' The actual per-species computation `collinearityCheckGerHabitat()` wraps
#' in `reproducible::Cache()` -- pulled into its own function so Cache()'s
#' digest covers exactly `spPa`/`blocksSp`/`requestedPredictors`/
#' `dropCollinearPredictors`/`spatialTerm` (what actually determines this
#' species' result), not the whole enclosing function's environment. A
#' change to just this species' predictor-table row (e.g. adding/removing
#' "hedges" for Neuntöter or Goldammer) changes `requestedPredictors`,
#' which changes the digest for exactly that species -- every other
#' species stays a cache hit.
#'
#' @param sp Character. Species Latin name.
#' @param spPa data.frame. This species' pooled occurrence+covariate data.
#' @param blocksSp List. This species' spatial-block-CV object.
#' @param requestedPredictors Character vector. This species'
#'   `speciesPredictorTable` entry -- the only source of its predictors.
#' @param dropCollinearPredictors Logical. Prune `requestedPredictors` for
#'   collinearity via `select07Blockcv()`, or use them exactly as given?
#' @param spatialTerm Logical. Append `x`/`y` as predictors?
#' @param threshold Numeric. Absolute correlation threshold.
#' @param univar Character. Initial univariate model form.
#' @return List with `data` (the final table) and `predictors`.
resolveHabitatSpeciesPredictors <- function(sp, spPa, blocksSp, requestedPredictors,
                                             dropCollinearPredictors, spatialTerm,
                                             threshold, univar) {
  requested <- requestedPredictors
  requested <- requested[!is.na(requested) & nzchar(trimws(requested))]

  if (length(requested) == 0) {
    stop(sprintf("Species '%s' has no predictors listed in speciesConfig_predictors.csv at HABITAT scale!", sp))
  }

  missingCols <- setdiff(requested, names(spPa))
  if (length(missingCols) > 0) {
    stop(sprintf("\n[PREDICTOR ERROR] Species '%s' [HABITAT scale]:\nRequested predictor(s) %s were NOT found in occurrence data!\nAvailable columns: %s\n",
                 sp, paste(dQuote(missingCols), collapse = ", "), paste(names(spPa), collapse = ", ")))
  }

  if (isTRUE(dropCollinearPredictors)) {
    nPres <- sum(spPa$occurrence == 1)
    nAbs <- sum(spPa$occurrence == 0)
    X <- spPa[, requested, drop = FALSE]
    varSel <- select07Blockcv(X = X, y = spPa$occurrence, threshold = threshold,
                              univar = univar, spBlock = blocksSp, weights = rep(1, nrow(spPa)))
    occNum <- max(floor(min(nPres, nAbs) / 10), 1)
    predSel <- stats::na.omit(varSel$pred_sel[1:min(occNum, length(varSel$pred_sel))])
  } else {
    predSel <- requested
  }

  # Append spatial coordinates ONCE if configured
  if (isTRUE(spatialTerm)) {
    predSel <- c(as.character(predSel), "x", "y")
  }

  message("Predictors used (", length(predSel), "): ", paste(predSel, collapse = ", "))

  spPaOut <- spPa
  spPaOut$foldID <- blocksSp$folds_ids

  keepCols <- c("AREA_NATIONAL_CODE", "year", "latin_name", "occurrence", "x", "y", "foldID", predSel)
  list(data = spPaOut[, intersect(keepCols, names(spPaOut))],
       predictors = as.character(predSel))
}
