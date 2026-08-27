#' Load per-species EBBA2 occurrence+bioclim data (European scale)
#'
#' European occurrence data has no year dimension to pool across -- one
#' RDS per species, as produced by dataPrep_Monitor's
#' `occurrencePrepEurope()`.
#'
#' @param occDir Character. Directory of `<species>_EBBA2_pa_env.rds` files.
#' @param species Character vector of Latin species names.
#' @return Named list (by species) of occurrence+bioclim data.frames.
poolOccurrenceEurope <- function(occDir, species) {
  occFiles <- list.files(occDir, pattern = "\\.rds$", full.names = TRUE)

  result <- list()
  for (sp in species) {
    spClean <- gsub(" ", "_", sp)
    spFile <- occFiles[grep(spClean, occFiles, fixed = TRUE)]
    if (length(spFile) == 0) {
      warning("No occurrence file found for: ", sp)
      next
    }
    result[[sp]] <- readRDS(spFile[1])
  }
  result
}
