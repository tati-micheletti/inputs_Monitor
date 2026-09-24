#' Candidate predictor column names for the German habitat and landscape SDMs
#'
#' Shared by both scales -- land use (13 categories) + land cover (3
#' categories) + DEM derivatives (3 layers).
#'
#' `hedges` deliberately excluded (methodology decision, 2026-09) -- still
#' computed upstream (dataPrep_Monitor's landuse category layers), just
#' never offered to collinearity selection, so no model ever uses it.
#'
#' @return Character vector of covariate column names.
covariatePredictorColumns <- function() {
  c("grassland", "winter_cereals", "summer_cereals", "maize", "root_crops",
    "rapeseed", "sunflower", "vegetables", "legumes", "fallow",
    "grapevine", "hops", "orchards_and_berries",
    "built_up", "trees", "water",
    "elevation", "slope", "solar_radiation")
}
