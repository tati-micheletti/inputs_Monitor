#' Non-spatial stand-in for spatial CV blocking
#'
#' Used in place of `spatialBlockingEurope()`/`spatialBlockingGerHabitat()`/
#' `spatialBlockingGerLandscape()` when `P(sim)$runSpatialBlocking` is
#' FALSE. Assigns a plain stratified random k-fold (stratified by
#' `occurrence` so every fold gets both presences and absences) instead
#' of a spatially-autocorrelation-aware block, but returns an object
#' with the SAME structural fields (`folds_ids`, `records`,
#' `blocks$block_id`) that `blockCV::cv_spatial()` returns -- so
#' `collinearityCheckX()` and any downstream code can treat a real or a
#' mimicked blocks object identically. This exists so you can compare
#' model performance with and without spatial (and, later, temporal)
#' blocking without the rest of the pipeline needing to change.
#'
#' @param paData data.frame with an `occurrence` (0/1) column.
#' @param k Integer. Number of folds, default 5.
#' @return List with `folds_ids`, `k`, `records`, `blocks` (block_id
#'   only), and `strategy = "mimic-random"`.
mimicSpatialBlocks <- function(paData, k = 5) {
  n <- nrow(paData)
  foldsIds <- integer(n)
  for (cls in unique(paData$occurrence)) {
    idx <- which(paData$occurrence == cls)
    foldsIds[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }

  records <- do.call(rbind, lapply(seq_len(k), function(fo) {
    data.frame(train_0 = sum(paData$occurrence[foldsIds != fo] == 0),
               train_1 = sum(paData$occurrence[foldsIds != fo] == 1),
               test_0 = sum(paData$occurrence[foldsIds == fo] == 0),
               test_1 = sum(paData$occurrence[foldsIds == fo] == 1))
  }))

  list(folds_ids = foldsIds,
       k = k,
       records = records,
       blocks = list(block_id = seq_len(k)),
       strategy = "mimic-random")
}
