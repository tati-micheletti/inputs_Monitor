#' Load and pool per-species-year German habitat occurrence data
#'
#' Same route appearing in multiple years is pooled into one table per
#' species (all years combined) -- a route from an unfiltered year still
#' correctly ends up in the same spatial block across years, avoiding
#' data leakage across CV folds.
#'
#' @param occDir Character. Directory of `<species>_habitat_<year>.rds` files.
#' @param species Character vector of Latin species names.
#' @param years Integer vector, or NULL (default). If supplied, only files
#'   whose embedded year is in this set are pooled -- confirmed necessary
#'   2026-09-30: without this, stale files left over from an earlier run
#'   with a different (e.g. broader) year range get silently pooled
#'   alongside the current run's fresh files, mixing different covariate
#'   vintages into one training table. NULL (default): pool every file
#'   found, today's original behavior.
#' @return Named list (by species) of pooled occurrence+covariate data.frames.
poolOccurrenceGerHabitat <- function(occDir, species, years = NULL) {
  occFiles <- list.files(occDir, pattern = "\\.rds$", full.names = TRUE)

  result <- list()
  for (sp in species) {
    spClean <- gsub(" ", "_", sp)
    spFiles <- occFiles[grep(spClean, occFiles, fixed = TRUE)]

    if (!is.null(years)) {
      fileYears <- as.integer(sub(".*_habitat_(\\d{4})\\.rds$", "\\1", basename(spFiles)))
      spFiles <- spFiles[fileYears %in% years]
    }

    if (length(spFiles) == 0) {
      message("No occurrence files found for ", sp, " -- skipping")
      next
    }
    spList <- lapply(spFiles, function(f) stripYearSuffixes(readRDS(f)))
    # dplyr::bind_rows() rather than rbind() -- see poolOccurrenceGerLandscape.R's
    # own docstring for the rationale (same pattern, same 2026-09-30 fix).
    result[[sp]] <- dplyr::bind_rows(spList)
  }
  result
}
