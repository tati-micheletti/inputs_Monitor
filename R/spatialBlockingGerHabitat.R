#' Create spatial CV blocks for the German habitat SDM (200m scale)
#'
#' Pooled approach: blocks are computed once per species on all years
#' combined, so the same route appearing in multiple years is always
#' assigned to the same block (avoids data leakage across folds). Block
#' size follows Wiedenroth et al.: `cv_spatial_autocor()` capped at
#' `maxBlockSizeM` (default 200km, empirically determined for the German
#' extent at 200m scale).
#'
#' @param pooledData Named list (by species) of pooled occurrence+covariate
#'   data.frames (see `poolOccurrenceGerHabitat()`).
#' @param refRasterPath Character. Path to a 200m reference raster (e.g.
#'   solar_radiation_habitat.tif) used as the `cv_spatial()` reference grid.
#' @param maxBlockSizeM Numeric. Maximum block size in metres.
#' @param k Integer. Number of folds, default 5.
#' @return Named list (by species) of `blockCV::cv_spatial()` result objects.
spatialBlockingGerHabitat <- function(pooledData, refRasterPath, maxBlockSizeM = 200000, k = 5) {

  refRaster <- terra::rast(refRasterPath)
  result <- list()

  for (sp in names(pooledData)) {
    message("\n  -- ", sp, " --------------------------")
    spPa <- pooledData[[sp]]

    nPres <- sum(spPa$occurrence == 1)
    nAbs <- sum(spPa$occurrence == 0)
    message("Pooled records: ", nrow(spPa), " (", nPres, " pres / ", nAbs, " abs)")
    message("Unique routes: ", length(unique(spPa$AREA_NATIONAL_CODE)))

    if (nPres < 10 || nAbs < 10) {
      message("Too few records -- skipping")
      next
    }

    sfOcc <- sf::st_as_sf(spPa, coords = c("x", "y"), crs = 3035)

    message("Determining block size...")
    cvBlSize <- determineBlockSize(sfOcc, maxBlockSizeM)

    message("  Creating spatial blocks (", k, " folds)...")
    scv <- tryCatch({
      blockCV::cv_spatial(x = sfOcc, column = "occurrence", r = refRaster, k = k,
                           size = cvBlSize, selection = "random", iteration = 50,
                           progress = FALSE, biomod2 = FALSE, plot = FALSE, report = FALSE)
    }, error = function(e) {
      warning("cv_spatial failed: ", e$message)
      NULL
    })
    if (is.null(scv)) next

    nBlocks <- length(scv$blocks$block_id)
    message("Blocks created: ", nBlocks)
    if (nBlocks < 15) {
      warning("Only ", nBlocks, " blocks for ", sp, " (minimum recommended: 15)")
    }

    zeroFolds <- which(scv$records$test_0 == 0 | scv$records$test_1 == 0)
    if (length(zeroFolds) > 0) {
      warning("Zero records in fold(s): ", paste(zeroFolds, collapse = ", "))
    }

    result[[sp]] <- scv
  }

  result
}
