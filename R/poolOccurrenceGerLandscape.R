#' Load and pool per-species-year German landscape occurrence data
#'
#' Same route appearing in multiple years is pooled into one table per
#' species (all years combined).
#'
#' @param occDir Character. Directory of `<species>_landscape_<year>.rds` files.
#' @param species Character vector of Latin species names.
#' @return Named list (by species) of pooled occurrence+covariate data.frames.
poolOccurrenceGerLandscape <- function(occDir, species) {
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
