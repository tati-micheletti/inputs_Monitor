#' Create spatial CV blocks for the German landscape SDM (1km scale)
#'
#' Pooled approach: blocks are computed once per species on all years
#' combined. Block size follows Wiedenroth et al.: `cv_spatial_autocor()`
#' capped at `maxBlockSizeM` (default 200km, empirically determined --
#' no indication in Wiedenroth et al. of what they used at this scale).
#'
#' @param pooledData Named list (by species) of pooled occurrence+covariate
#'   data.frames (see `poolOccurrenceGerLandscape()`).
#' @param refRasterPath Character. Path to a 1km reference raster (e.g.
#'   the first band of a landuse_<year>_landscape.tif) used as the
#'   `cv_spatial()` reference grid.
#' @param maxBlockSizeM Numeric. Maximum block size in metres.
#' @param k Integer. Number of folds, default 5.
#' @return Named list (by species) of `blockCV::cv_spatial()` result objects.
spatialBlockingGerLandscape <- function(pooledData, refRasterPath, maxBlockSizeM = 200000, k = 5) {

  refRaster <- terra::rast(refRasterPath)[[1]]
  result <- list()

  for (sp in names(pooledData)) {
    message("\n  -- ", sp, " --------------------")
    spPa <- pooledData[[sp]]

    message("Pooled records: ", nrow(spPa), " (", sum(spPa$occurrence == 1), " pres / ",
            sum(spPa$occurrence == 0), " abs)")
    message("Unique routes: ", length(unique(spPa$ROUTENCODE)))

    nPres <- sum(spPa$occurrence == 1)
    nAbs <- sum(spPa$occurrence == 0)
    if (nPres < 10 || nAbs < 10) {
      message("Too few records -- skipping")
      next
    }

    sfOcc <- sf::st_as_sf(spPa, coords = c("x", "y"), crs = 3035)

    message("  Determining block size...")
    cvBlSize <- determineBlockSize(sfOcc, maxBlockSizeM)

    message("  Creating spatial blocks (", k, " folds)...")
    scv <- tryCatch({
      blockCV::cv_spatial(x = sfOcc, column = "occurrence", r = refRaster, k = k,
                           size = cvBlSize, selection = "random", iteration = 50,
                           progress = FALSE, biomod2 = FALSE, plot = FALSE, report = FALSE)
    }, error = function(e) {
      warning("  cv_spatial failed: ", e$message)
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
      warning("  Zero records in fold(s): ", paste(zeroFolds, collapse = ", "))
    }

    result[[sp]] <- scv
  }

  result
}
