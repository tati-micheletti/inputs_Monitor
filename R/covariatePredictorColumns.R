#' Candidate predictor column names for the German habitat and landscape SDMs
#'
#' Shared by both scales -- land use (14 categories) + land cover (3
#' categories) + DEM derivatives (3 layers).
#'
#' @return Character vector of covariate column names.
covariatePredictorColumns <- function() {
  c("grassland", "winter_cereals", "summer_cereals", "maize", "root_crops",
    "rapeseed", "sunflower", "vegetables", "legumes", "hedges", "fallow",
    "grapevine", "hops", "orchards_and_berries",
    "built_up", "trees", "water",
    "elevation", "slope", "solar_radiation")
}
