#' Load and pool per-species-year German habitat occurrence data
#'
#' Same route appearing in multiple years is pooled into one table per
#' species (all years combined) -- a route from an unfiltered year still
#' correctly ends up in the same spatial block across years, avoiding
#' data leakage across CV folds.
#'
#' @param occDir Character. Directory of `<species>_habitat_<year>.rds` files.
#' @param species Character vector of Latin species names.
#' @return Named list (by species) of pooled occurrence+covariate data.frames.
poolOccurrenceGerHabitat <- function(occDir, species) {
  occFiles <- list.files(occDir, pattern = "\\.rds$", full.names = TRUE)

  result <- list()
  for (sp in species) {
    spClean <- gsub(" ", "_", sp)
    spFiles <- occFiles[grep(spClean, occFiles, fixed = TRUE)]
    if (length(spFiles) == 0) {
      message("No occurrence files found for ", sp, " -- skipping")
      next
    }
    spList <- lapply(spFiles, function(f) stripYearSuffixes(readRDS(f)))
    result[[sp]] <- do.call(rbind, spList)
  }
  result
}
