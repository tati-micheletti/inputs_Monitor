#' Create spatial CV blocks for the European climate SDM
#'
#' Block size follows Wiedenroth et al.: `cv_spatial_autocor()` for the
#' natural spatial-autocorrelation range, capped at `maxBlockSizeM`
#' (default 1500km, tuned for >=15 blocks across the European EBBA2
#' extent) and floored at `minBlockSizeM` (default 200km, 2x the
#' European thinning distance).
#'
#' @param pooledData Named list (by species) of occurrence+bioclim data.frames.
#' @param bioclimFile Character. Path to the training-year bioclim raster
#'   (used as the `cv_spatial()` reference grid).
#' @param maxBlockSizeM Numeric. Maximum block size in metres.
#' @param minBlockSizeM Numeric. Minimum block size in metres.
#' @param k Integer. Number of folds, default 5.
#' @return Named list (by species) of `blockCV::cv_spatial()` result objects.
spatialBlockingEurope <- function(pooledData, bioclimFile, maxBlockSizeM = 1500000,
                                   minBlockSizeM = 200000, k = 5) {

  bioclim <- terra::rast(bioclimFile)
  result <- list()

  for (sp in names(pooledData)) {
    message("\n  Processing: ", sp)
    spPa <- pooledData[[sp]]
    message("    Records: ", nrow(spPa), " (", sum(spPa$occurrence == 1), " pres / ",
            sum(spPa$occurrence == 0), " abs)")

    sfSpecOcc <- sf::st_as_sf(spPa, coords = c("x", "y"), crs = terra::crs(bioclim))

    message("    Determining block size via spatial autocorrelation...")
    cvBlSize <- determineBlockSize(sfSpecOcc, maxBlockSizeM, minBlockSizeM)

    message("    Creating spatial blocks (", k, " folds)...")
    scv <- blockCV::cv_spatial(x = sfSpecOcc, column = "occurrence", r = bioclim, k = k,
                                size = cvBlSize, selection = "random", iteration = 50,
                                progress = FALSE, biomod2 = FALSE, plot = FALSE, report = FALSE)

    nBlocks <- length(scv$blocks$block_id)
    message("    Blocks created: ", nBlocks)
    if (nBlocks < 15) {
      warning("    Only ", nBlocks, " blocks created for ", sp,
              " (minimum recommended is 15). Consider reducing block size manually.")
    }

    foldRecords <- scv$records
    if (!is.null(foldRecords)) {
      for (foldI in seq_len(nrow(foldRecords))) {
        testPres <- foldRecords[foldI, "test_0"]
        testAbs <- foldRecords[foldI, "test_1"]
        if (any(c(testPres, testAbs) == 0)) {
          message("    Fold ", foldI, " has zero records for a class -- presences in test: ",
                  testPres, ", absences in test: ", testAbs)
        }
      }
    }

    result[[sp]] <- scv
  }

  result
}
