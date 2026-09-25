#' Pull one scale's slice out of a per-species, per-scale nested structure
#'
#' `speciesConfig_predictors.csv`'s loader (`loadSpeciesPredictorConfig()` in
#' sharedSpeciesConfig.R) produces a structure keyed species -> scale ->
#' character vector (which predictor names that species includes at that
#' scale). Each scale's `collinearityCheck*()` function only needs its own
#' scale's slice, keyed just by species, to pass as its `predictorsToUse`
#' argument -- this reshapes that once per scale instead of each function
#' reaching into the nested structure itself.
#'
#' @param perSpeciesNested The nested list, or NULL.
#' @param scale Character, one of "climate"/"landscape"/"habitat".
#' @return Named list (species -> that scale's character vector), or NULL if
#'   `perSpeciesNested` is NULL. A species with nothing for this scale is
#'   simply absent from the result (falls through to that function's normal
#'   `runCollinearityCheck` default, same as any other unlisted species).
extractScaleExtras <- function(perSpeciesNested, scale) {
  if (is.null(perSpeciesNested)) return(NULL)
  out <- lapply(perSpeciesNested, function(sp) sp[[scale]])
  out[!sapply(out, is.null)]
}
