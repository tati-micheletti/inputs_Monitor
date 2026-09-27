#' Candidate predictor column names for the German habitat and landscape SDMs
#'
#' Shared by both scales -- land use (14 categories) + land cover 
#' (3 categories) + DEM derivatives (2
#' layers: elevation, slope -- solar_radiation dropped 2026-09-26, see
#' DECISIONS.md).
#'
#' `hedges` is included by default (methodology decision, 2026-09), but can be
#' excluded by marking in the shared config files. Either way it's already
#' computed upstream in dataPrep_Monitor's landuse category layers, with
#' pre-2017/2022-2023 gaps already backfilled from the nearest real year
#' (`loadCovariates()`/`loadHabitatCovariates()`/`occurrencePrepGerHabitat()`
#' in dataPrep_Monitor) -- "backfill" here doesn't invent new fill logic, it
#' just re-offers an already-backfilled column to collinearity selection.
#' `collinearityCheckGerHabitat()`/`GerLandscape()` still drop it per-species
#' if it ends up all-NA for that species' pooled data, same as any other
#' candidate predictor.
#'
#' @return Character vector of covariate column names.
covariatePredictorColumns <- function(...) {
  c("grassland", "winter_cereals", "summer_cereals", "maize", "root_crops",
    "rapeseed", "sunflower", "vegetables", "legumes", "hedges", "fallow",
    "grapevine", "hops", "orchards_and_berries",
    "built_up", "trees", "water",
    "elevation", "slope")
}
