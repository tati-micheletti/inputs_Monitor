#' Candidate predictor column names for the German habitat and landscape SDMs
#'
#' Shared by both scales -- land use (14 categories) + land cover
#' (7 categories) + DEM derivatives (2
#' layers: elevation, slope -- solar_radiation dropped 2026-09-26, see
#' DECISIONS.md) + two derived covariates (dist_to_woodland,
#' landscape_heterogeneity -- see dataPrep_Monitor's computeDistToWoodland()/
#' computeHeterogeneityIndex(), added 2026-09-29 as hedges-backfill
#' alternatives).
#'
#' `hedges` is included by default (methodology decision, 2026-09), but can be
#' excluded by marking in the shared config files. Either way it's already
#' computed upstream in dataPrep_Monitor's landuse category layers, and since
#' 2026-11-03 it is ONE FIXED layer in every land-use file (the mean of the hedge
#' layers of the years that have CTM hedge data, 2017-2021 and 2024-2025;
#' `completeHedges(mode = "mean")`, see DECISIONS.md 2026-11-03), so hedges carry
#' no time signal. (Before that, years without hedge data took the layer of the
#' nearest preceding year; that rule is kept as `mode = "nearest"` and survives,
#' unused, as back-fill code in `loadCovariates()`/`loadHabitatCovariates()`
#' and `occurrencePrepGerHabitat()`, which no longer triggers because the files
#' have no NA hedge cells.) "Backfill" in the old text of this function only
#' re-offered an already-filled column to collinearity selection.
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
    "natural_grassland", "heath_moor", "scrub", "wetland_peat",         # CORINE semi-natural classes (2026-10-21)
    "elevation", "slope",
    "dist_to_woodland", "landscape_heterogeneity", "unmapped_other")
}
