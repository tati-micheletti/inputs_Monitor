#' Resolve final predictors and assemble the model-ready table (Europe scale)
#'
#' For every species, always saves a correlation-structure corrplot over
#' all candidate bioclim variables (diagnostic, independent of the
#' predictor-resolution strategy below), then resolves the predictor set
#' for the final table using, in priority order:
#'   1. `predictorsToUse` if not NULL (an explicit vector, or "all") --
#'      always wins, with a warning if `runCollinearityCheck` was TRUE.
#'   2. Real block-CV collinearity selection (`select07Blockcv()`) if
#'      `runCollinearityCheck` is TRUE, capped at 1 predictor per 10
#'      occurrences (the rarer of presence/absence counts).
#'   3. All available bioclim columns, unfiltered, otherwise.
#'
#' Regardless of which branch resolves the predictors, the output table
#' always has the same guaranteed base columns (cell50x50,
#' birdlife_scientific_name, occurrence, x, y, foldID) plus whatever
#' predictor columns were resolved -- this is what keeps
#' `runSpatialBlocking`/`runCollinearityCheck` safely toggleable without
#' breaking models_Monitor downstream.
#'
#' @param pooledData Named list (by species) of occurrence+bioclim data.frames.
#' @param blocksData Named list (by species) of blocks objects (real or
#'   mimicked) -- must have `$folds_ids` aligned to `pooledData[[sp]]` rows.
#' @param runCollinearityCheck Logical. Whether to run real collinearity selection.
#' @param predictorsToUse NULL, "all", or a character vector of predictor names.
#' @param corrplotDir Character. Directory to save correlation plots in.
#' @param threshold Numeric. Absolute correlation threshold, default 0.7.
#' @param univar Character. Initial univariate model form, default "gam".
#' @param perSpeciesExtraCandidates Named list, or NULL (default). Per-species
#'   extra candidate predictor names (character vectors) to ADD to the
#'   default bioclim candidate pool before collinearity selection runs --
#'   see `loadSpeciesPredictorExtras()`. Does not replace the default pool.
#' @return Named list (by species) with `data` (the final table) and
#'   `predictors` (character vector of predictor columns used).
collinearityCheckEurope <- function(pooledData, blocksData, runCollinearityCheck,
                                     predictorsToUse, corrplotDir, threshold = 0.7,
                                     univar = "gam", perSpeciesExtraCandidates = NULL) {

  dir.create(corrplotDir, recursive = TRUE, showWarnings = FALSE)
  result <- list()

  for (sp in names(pooledData)) {
    message("\n  Processing: ", sp)
    spPa <- pooledData[[sp]]
    spClean <- gsub(" ", "_", sp)

    bioVars <- bioclimPredictorColumns()
    missing <- setdiff(bioVars, names(spPa))
    if (length(missing) > 0) {
      stop("Missing bioclim variables for ", sp, ": ", paste(missing, collapse = ", "))
    }

    if (!is.null(perSpeciesExtraCandidates) && sp %in% names(perSpeciesExtraCandidates)) {
      extras <- perSpeciesExtraCandidates[[sp]]
      missingExtras <- setdiff(extras, names(spPa))
      if (length(missingExtras) > 0) {
        warning(sp, ": extra candidate predictor(s) not found in the data, ignored: ",
                paste(missingExtras, collapse = ", "))
      }
      bioVars <- union(bioVars, intersect(extras, names(spPa)))
    }

    corMat <- cor(spPa[, bioVars], method = "spearman")
    corrplotFile <- file.path(corrplotDir, paste0(spClean, "_EU.png"))
    grDevices::png(corrplotFile, width = 1600, height = 1600, res = 240)
    corrplot::corrplot.mixed(corMat, tl.pos = "lt", tl.cex = 0.6, number.cex = 0.5,
                              addCoefasPercent = TRUE)
    grDevices::dev.off()

    nPres <- sum(spPa$occurrence == 1)
    nAbs <- sum(spPa$occurrence == 0)

    if (!is.null(predictorsToUse)) {
      if (isTRUE(runCollinearityCheck)) {
        warning("predictorsToUse overrides runCollinearityCheck=TRUE for ", sp,
                " -- using the specified predictor set instead of collinearity selection.")
      }
      predSel <- if (identical(predictorsToUse, "all")) bioVars else intersect(predictorsToUse, bioVars)
    } else if (isTRUE(runCollinearityCheck)) {
      blocksSp <- blocksData[[sp]]
      varSel <- select07Blockcv(X = spPa[, bioVars], y = spPa$occurrence, threshold = threshold,
                                 univar = univar, spBlock = blocksSp, weights = rep(1, nrow(spPa)))
      occNum <- max(floor(min(nPres, nAbs) / 10), 1)
      predSel <- stats::na.omit(varSel$pred_sel[1:min(occNum, length(varSel$pred_sel))])
    } else {
      predSel <- bioVars
    }

    message("Predictors used (", length(predSel), "): ", paste(predSel, collapse = ", "))

    spPaOut <- spPa
    spPaOut$foldID <- blocksData[[sp]]$folds_ids

    keepCols <- c("cell50x50", "birdlife_scientific_name", "occurrence", "x", "y", "foldID", predSel)
    result[[sp]] <- list(data = spPaOut[, intersect(keepCols, names(spPaOut))],
                          predictors = as.character(predSel))
  }

  result
}
