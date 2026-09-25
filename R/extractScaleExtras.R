#' Pull one scale's slice out of a per-species, per-scale nested structure
#'
#' Both `perSpeciesHedges` and `perSpeciesExtraCandidates` (see
#' `inputs_Monitor`'s parameters of the same names, sourced from
#' `loadSpeciesGeneralConfig()`/`loadSpeciesPredictorExtras()`) are keyed
#' species -> scale -> value (a "drop"/"backfill" string for the former, a
#' character vector of extra predictor names for the latter). Each scale's
#' `collinearityCheck*()` function only needs its own scale's slice, keyed
#' just by species -- this reshapes that once per scale instead of each
#' function reaching into the nested structure itself. Generic over what
#' the per-species value actually is (string or vector), so it serves both.
#'
#' @param perSpeciesNested The nested list, or NULL.
#' @param scale Character, one of "climate"/"landscape"/"habitat".
#' @return Named list (species -> that scale's value), or NULL if
#'   `perSpeciesNested` is NULL. A species with nothing for this scale is
#'   simply absent from the result.
extractScaleExtras <- function(perSpeciesNested, scale) {
  if (is.null(perSpeciesNested)) return(NULL)
  out <- lapply(perSpeciesNested, function(sp) sp[[scale]])
  out[!sapply(out, is.null)]
}
