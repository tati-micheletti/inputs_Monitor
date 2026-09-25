#' Resolve final predictors and assemble the model-ready table (German landscape scale)
#'
#' Same predictor-resolution priority as `collinearityCheckEurope()`:
#' `predictorsToUse` override > real collinearity selection > all
#' available covariates -- guaranteeing the same base output columns
#' (ROUTENCODE, latin_name, occurrence, x, y, foldID) regardless of the
#' `runSpatialBlocking`/`runCollinearityCheck` toggles.
#'
#' @param pooledData Named list (by species) of pooled occurrence+covariate
#'   data.frames (see `poolOccurrenceGerLandscape()`).
#' @param blocksData Named list (by species) of blocks objects (real or
#'   mimicked) -- must have `$folds_ids` aligned to `pooledData[[sp]]` rows.
#' @param runCollinearityCheck Logical. Whether to run real collinearity selection.
#' @param predictorsToUse NULL, "all", or a character vector of predictor names.
#' @param corrplotDir Character. Directory to save correlation plots in.
#' @param threshold Numeric. Absolute correlation threshold, default 0.7.
#' @param univar Character. Initial univariate model form, default "gam".
#' @param hedgesTreatment Character, "drop" (default) or "backfill" -- see
#'   `covariatePredictorColumns()`. Applies to every species unless
#'   overridden per-species in `perSpeciesHedges`.
#' @param perSpeciesHedges Named character vector/list, or NULL (default).
#'   Per-species `hedgesTreatment` overrides, keyed by species Latin name.
#' @param perSpeciesExtraCandidates Named list, or NULL (default). Per-species
#'   extra candidate predictor names (character vectors) to ADD to this
#'   scale's default candidate pool before collinearity selection runs --
#'   see `loadSpeciesPredictorExtras()`. Does not replace the default pool.
#' @return Named list (by species) with `data` (the final table) and
#'   `predictors` (character vector of predictor columns used).
collinearityCheckGerLandscape <- function(pooledData, blocksData, runCollinearityCheck,
                                           predictorsToUse, corrplotDir, threshold = 0.7,
                                           univar = "gam", hedgesTreatment = "drop",
                                           perSpeciesHedges = NULL, perSpeciesExtraCandidates = NULL) {

  dir.create(corrplotDir, recursive = TRUE, showWarnings = FALSE)
  result <- list()

  for (sp in names(pooledData)) {
    message("\n  -- ", sp, " --------------------")
    spPa <- pooledData[[sp]]
    spClean <- gsub(" ", "_", sp)

    if (is.null(blocksData[[sp]])) {
      message("No blocks available for ", sp, " -- skipping")
      next
    }

    spHedgesTreatment <- if (!is.null(perSpeciesHedges) && sp %in% names(perSpeciesHedges)) {
      perSpeciesHedges[[sp]]
    } else {
      hedgesTreatment
    }
    allPredictors <- covariatePredictorColumns(spHedgesTreatment)
    if (!is.null(perSpeciesExtraCandidates) && sp %in% names(perSpeciesExtraCandidates)) {
      allPredictors <- union(allPredictors, perSpeciesExtraCandidates[[sp]])
    }

    nPres <- sum(spPa$occurrence == 1)
    nAbs <- sum(spPa$occurrence == 0)
    message("Pooled records: ", nrow(spPa), " (", nPres, " pres / ", nAbs, " abs)")

    predCols <- intersect(allPredictors, names(spPa))
    allNaCols <- predCols[sapply(predCols, function(col) all(is.na(spPa[[col]])))]
    if (length(allNaCols) > 0) {
      message("Removing all-NA predictors: ", paste(allNaCols, collapse = ", "))
      predCols <- setdiff(predCols, allNaCols)
    }

    X <- spPa[, predCols, drop = FALSE]
    X <- X[, sapply(X, function(col) length(unique(col[!is.na(col)])) > 1), drop = FALSE]

    if (ncol(X) < 2) {
      warning("Too few variable predictors for ", sp, " -- skipping")
      next
    }

    corrplotFile <- file.path(corrplotDir, paste0(spClean, "_landscape_pooled.png"))
    if (!file.exists(corrplotFile)) {
      grDevices::png(corrplotFile, width = 1200, height = 1200, res = 150)
      corrplot::corrplot.mixed(cor(X, method = "spearman", use = "complete.obs"),
                                tl.pos = "lt", tl.cex = 0.6, number.cex = 0.4, addCoefasPercent = TRUE)
      grDevices::dev.off()
    }

    if (!is.null(predictorsToUse)) {
      if (isTRUE(runCollinearityCheck)) {
        warning("predictorsToUse overrides runCollinearityCheck=TRUE for ", sp,
                " -- using the specified predictor set instead of collinearity selection.")
      }
      predSel <- if (identical(predictorsToUse, "all")) colnames(X) else intersect(predictorsToUse, colnames(X))
    } else if (isTRUE(runCollinearityCheck)) {
      varSel <- tryCatch({
        select07Blockcv(X = X, y = spPa$occurrence, threshold = threshold, univar = univar,
                         spBlock = blocksData[[sp]], weights = rep(1, nrow(spPa)))
      }, error = function(e) {
        warning("Variable selection failed for ", sp, ": ", e$message)
        NULL
      })
      if (is.null(varSel)) next
      occNum <- max(floor(min(nPres, nAbs) / 10), 1)
      predSel <- stats::na.omit(varSel$pred_sel[1:min(occNum, length(varSel$pred_sel))])
    } else {
      predSel <- colnames(X)
    }

    message("Predictors used (", length(predSel), "): ", paste(predSel, collapse = ", "))

    spPaOut <- spPa
    spPaOut$foldID <- blocksData[[sp]]$folds_ids

    keepCols <- c("ROUTENCODE", "latin_name", "occurrence", "x", "y", "foldID", predSel)
    result[[sp]] <- list(data = spPaOut[, intersect(keepCols, names(spPaOut))],
                          predictors = as.character(predSel))
  }

  result
}
