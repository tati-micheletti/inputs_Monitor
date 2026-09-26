#' Resolve final predictors and assemble the model-ready table (German habitat scale)
#'
#' Three predictor-resolution modes, chosen via `predictorsToUse` (per
#' species, or one value for everyone): `"table"` (use
#' `speciesPredictorTable`'s exact list for that species), `"all"` (use
#' every available covariate, unfiltered), or `"auto"` (real block-CV
#' collinearity selection via `select07Blockcv()`, capped at 1 predictor
#' per 10 occurrences). This guarantees the same base output columns
#' (AREA_NATIONAL_CODE, year, latin_name, occurrence, x, y, foldID)
#' regardless of which mode a species uses.
#'
#' @param pooledData Named list (by species) of pooled occurrence+covariate
#'   data.frames (see `poolOccurrenceGerHabitat()`).
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
#' @param hedgesTreatment Character, "drop" (default) or "backfill" -- see
#'   `covariatePredictorColumns()`. A single shared value, tuned in code, not
#'   per-species -- per-species hedges inclusion belongs in
#'   `speciesPredictorTable` instead (list "hedges" for whichever species
#'   should get it, once this is set to "backfill" so the column actually
#'   exists to list).
#' @return Named list (by species) with `data` (the final table) and
#'   `predictors` (character vector of predictor columns used -- always
#'   includes `x`/`y` on top of whatever the mode resolves, see
#'   DECISIONS.md's 2026-09-26 spatial-coordinate-predictor entry).
collinearityCheckGerHabitat <- function(pooledData, blocksData, predictorsToUse,
                                         speciesPredictorTable = NULL, corrplotDir,
                                         threshold = 0.7, univar = "gam",
                                         hedgesTreatment = "drop") {

  dir.create(corrplotDir, recursive = TRUE, showWarnings = FALSE)
  allPredictors <- covariatePredictorColumns(hedgesTreatment)
  result <- list()

  for (sp in names(pooledData)) {
    message("\n  -- ", sp, " --------------------------")
    spPa <- pooledData[[sp]]
    spClean <- gsub(" ", "_", sp)

    if (is.null(blocksData[[sp]])) {
      message("No blocks available for ", sp, " -- skipping")
      next
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

    corrplotFile <- file.path(corrplotDir, paste0(spClean, "_habitat_pooled.png"))
    if (!file.exists(corrplotFile)) {
      grDevices::png(corrplotFile, width = 1200, height = 1200, res = 150)
      corrplot::corrplot.mixed(cor(X, method = "spearman", use = "complete.obs"),
                                tl.pos = "lt", tl.cex = 0.6, number.cex = 0.4, addCoefasPercent = TRUE)
      grDevices::dev.off()
    }

    if (identical(spMode, "table") && is.null(speciesPredictorTable[[sp]])) {
      warning(sp, ": predictor_mode is \"table\" but no entry exists in ",
              "speciesPredictorTable -- falling back to \"auto\".")
      spMode <- "auto"
    }

    if (identical(spMode, "table")) {
      predSel <- intersect(speciesPredictorTable[[sp]], colnames(X))
    } else if (identical(spMode, "all")) {
      predSel <- colnames(X)
    } else {
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

    keepCols <- c("AREA_NATIONAL_CODE", "year", "latin_name", "occurrence", "x", "y", "foldID", predSel)
    result[[sp]] <- list(data = spPaOut[, intersect(keepCols, names(spPaOut))],
                          predictors = as.character(predSel))
  }

  result
}
