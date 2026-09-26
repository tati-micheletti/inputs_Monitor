#' Resolve final predictors and assemble the model-ready table (Europe scale)
#'
#' For every species, always saves a correlation-structure corrplot over
#' all candidate bioclim variables (diagnostic, independent of the
#' predictor-resolution mode below). Three predictor-resolution modes,
#' chosen via `predictorsToUse` (per species, or one value for everyone):
#' `"table"` (use `speciesPredictorTable`'s exact list for that species),
#' `"all"` (use every bioclim variable, unfiltered), or `"auto"` (real
#' block-CV collinearity selection via `select07Blockcv()`, capped at 1
#' predictor per 10 occurrences).
#'
#' Regardless of which mode resolves the predictors, the output table
#' always has the same guaranteed base columns (cell50x50,
#' birdlife_scientific_name, occurrence, x, y, foldID) plus whatever
#' predictor columns were resolved.
#'
#' @param pooledData Named list (by species) of occurrence+bioclim data.frames.
#' @param blocksData Named list (by species) of blocks objects (real or
#'   mimicked) -- must have `$folds_ids` aligned to `pooledData[[sp]]` rows.
#' @param predictorsToUse Character `"table"`/`"all"`/`"auto"` (applied to
#'   every species), OR a named list (species -> one of those 3 strings) for
#'   per-species modes -- e.g. sourced from `speciesConfig_general.csv`'s
#'   `predictor_mode` column. A species absent from the list defaults to
#'   `"auto"`.
#' @param speciesPredictorTable Named list (species -> character vector), or
#'   NULL. Only consulted for species in `"table"` mode -- e.g. sourced from
#'   `speciesConfig_predictors.csv` via `loadSpeciesPredictorConfig()`. A
#'   `"table"`-mode species missing here falls back to `"auto"` with a
#'   warning, rather than silently using nothing.
#' @param corrplotDir Character. Directory to save correlation plots in.
#' @param threshold Numeric. Absolute correlation threshold, default 0.7.
#' @param univar Character. Initial univariate model form, default "gam".
#' @return Named list (by species) with `data` (the final table) and
#'   `predictors` (character vector of predictor columns used -- always
#'   includes `x`/`y` on top of whatever the mode resolves, see
#'   DECISIONS.md's 2026-09-26 spatial-coordinate-predictor entry).
collinearityCheckEurope <- function(pooledData, blocksData, predictorsToUse,
                                     speciesPredictorTable = NULL, corrplotDir,
                                     threshold = 0.7, univar = "gam") {

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

    spMode <- if (is.list(predictorsToUse)) {
      if (sp %in% names(predictorsToUse)) predictorsToUse[[sp]] else "auto"
    } else {
      predictorsToUse
    }
    if (!spMode %in% c("table", "all", "auto")) {
      stop(sp, ": invalid predictorsToUse mode '", spMode, "' -- must be ",
           "\"table\", \"all\", or \"auto\".")
    }

    corMat <- cor(spPa[, bioVars], method = "spearman")
    corrplotFile <- file.path(corrplotDir, paste0(spClean, "_EU.png"))
    grDevices::png(corrplotFile, width = 1600, height = 1600, res = 240)
    corrplot::corrplot.mixed(corMat, tl.pos = "lt", tl.cex = 0.6, number.cex = 0.5,
                              addCoefasPercent = TRUE)
    grDevices::dev.off()

    nPres <- sum(spPa$occurrence == 1)
    nAbs <- sum(spPa$occurrence == 0)

    if (identical(spMode, "table") && is.null(speciesPredictorTable[[sp]])) {
      warning(sp, ": predictor_mode is \"table\" but no entry exists in ",
              "speciesPredictorTable -- falling back to \"auto\".")
      spMode <- "auto"
    }

    if (identical(spMode, "table")) {
      predSel <- intersect(speciesPredictorTable[[sp]], bioVars)
    } else if (identical(spMode, "all")) {
      predSel <- bioVars
    } else {
      blocksSp <- blocksData[[sp]]
      varSel <- select07Blockcv(X = spPa[, bioVars], y = spPa$occurrence, threshold = threshold,
                                 univar = univar, spBlock = blocksSp, weights = rep(1, nrow(spPa)))
      occNum <- max(floor(min(nPres, nAbs) / 10), 1)
      predSel <- stats::na.omit(varSel$pred_sel[1:min(occNum, length(varSel$pred_sel))])
    }

    # Always add projected x/y coordinates as predictors, regardless of
    # mode -- a spatial trend-surface term meant to soak up residual
    # regional structure the environmental covariates alone don't capture
    # (see DECISIONS.md, 2026-09-26). Not a formal random effect (BRT/
    # dismo::gbm.step() has no mixed-model machinery) -- functionally,
    # letting the tree split on location itself.
    predSel <- c(as.character(predSel), "x", "y")

    message("Predictors used (", length(predSel), "): ", paste(predSel, collapse = ", "))

    spPaOut <- spPa
    spPaOut$foldID <- blocksData[[sp]]$folds_ids

    keepCols <- c("cell50x50", "birdlife_scientific_name", "occurrence", "x", "y", "foldID", predSel)
    result[[sp]] <- list(data = spPaOut[, intersect(keepCols, names(spPaOut))],
                          predictors = as.character(predSel))
  }

  result
}
