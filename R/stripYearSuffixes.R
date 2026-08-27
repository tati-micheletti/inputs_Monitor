#' Strip trailing `_<4-digit-year>` suffixes from column names
#'
#' Defensive normalization applied when pooling per-year occurrence
#' files across years, so that a covariate column named consistently
#' (e.g. "grassland") is not accidentally duplicated as "grassland_2022",
#' "grassland_2023", etc. across pooled years.
#'
#' @param df data.frame.
#' @return `df` with column names stripped of any trailing year suffix.
stripYearSuffixes <- function(df) {
  names(df) <- gsub("_\\d{4}$", "", names(df))
  df
}
