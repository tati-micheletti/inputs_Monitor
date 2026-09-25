defineModule(sim, list(
  name = "inputs_Monitor",
  description = paste("Creates the final per-species analysis tables for the bird monitor",
                       "pipeline: spatial CV blocking (real or a non-spatial mimic) and",
                       "collinearity-based predictor selection (or a user override), producing",
                       "the model-ready input tables consumed by models_Monitor."),
  keywords = c("bird monitor", "spatial blocking", "collinearity", "predictor selection"),
  authors = structure(list(list(given = "Tati", family = "Micheletti", role = c("aut", "cre"),
                                 email = "tati.micheletti@gmail.com", comment = NULL),
                           list(given = "Lisa", family = "Hildebrand", role = "aut",
                                email = "lisa.hildebrand@ufz.de", comment = NULL)), class = "person"),
  childModules = character(0),
  version = list(inputs_Monitor = "0.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("NEWS.md", "README.md", "inputs_Monitor.Rmd"),
  reqdPkgs = list("PredictiveEcology/SpaDES.core@development (>= 3.2.0)",
                   "terra", "sf", "blockCV", "dismo", "mgcv", "corrplot"),
  parameters = bindrows(
    defineParameter(".plots", "character", "screen", NA, NA,
                    "Used by Plots function, which can be optionally used here"),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA,
                    "Describes the simulation time at which the first plot event should occur."),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "Describes the simulation time interval between plot events."),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "Describes the simulation time at which the first save event should occur."),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events."),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used - e.g., a hash of the study",
                          "area obtained using `reproducible::studyAreaName()`"),
    ## .seed is optional: `list('init' = 123)` will `set.seed(123)` for the `init` event only.
    defineParameter(".seed", "list", list(), NA, NA,
                    "Named list of seeds to use for each event (names)."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Should caching of events or module be used?"),

    ## Toggles ---------------------------------------------------------------------
    ## NOTE: both spatialBlocking and collinearityCheck ALWAYS run every simulation
    ## -- these parameters control which STRATEGY they use internally (real vs a
    ## structurally-identical fallback), not whether they run at all. This is what
    ## guarantees sim$inputsData always has the same base table structure regardless
    ## of how these are set: see mimicSpatialBlocks() and collinearityCheckEurope()
    ## (and its ger_habitat/ger_landscape siblings) for exactly how each fallback
    ## keeps the contract.
    defineParameter("runSpatialBlocking", "logical", TRUE, NA, NA,
                    "TRUE: use real blockCV::cv_spatial() spatial CV blocks (Wiedenroth et al.",
                    "methodology). FALSE: use mimicSpatialBlocks() instead -- a plain stratified",
                    "random k-fold with the exact same output structure, so you can compare model",
                    "performance with/without spatial-autocorrelation-aware blocking."),
    defineParameter("predictorsToUse", "character", "auto", NA, NA,
                    "One of \"table\" (use speciesPredictorTable's exact list), \"all\" (use every",
                    "available covariate, unfiltered), or \"auto\" (default -- real block-CV",
                    "collinearity selection via select07Blockcv(), Wiedenroth et al. methodology,",
                    "capped at 1 predictor per 10 occurrences). Applied to every species, UNLESS",
                    "overridden per-species by a named list (species -> one of those 3 strings) --",
                    "e.g. sourced from speciesConfig_general.csv's predictor_mode column via",
                    "loadSpeciesGeneralConfig() in sharedSpeciesConfig.R, resolved once by",
                    "runMe.R/the orchestrating script and passed in here as a plain value. A",
                    "species absent from that list defaults to \"auto\"."),
    defineParameter("speciesPredictorTable", "list", NULL, NA, NA,
                    "NULL (default): no species can use predictorsToUse=\"table\" mode (falls back",
                    "to \"auto\" with a warning if any does). Otherwise a named list, species ->",
                    "scale -> character vector, the exact predictor set a \"table\"-mode species",
                    "uses at that scale. Sourced from speciesConfig_predictors.csv (repo root) via",
                    "loadSpeciesPredictorConfig() in sharedSpeciesConfig.R -- resolved once by the",
                    "orchestrating script and passed in as a plain value, same pattern as",
                    "sharedConfig.R's other shared values; this module doesn't read the CSV",
                    "itself, to stay self-contained/portable off this repo layout."),
    defineParameter("hedgesTreatment", "character", "drop", NA, NA,
                    "One of \"drop\" (default) or \"backfill\". \"drop\": hedges is never offered",
                    "to collinearity selection (methodology decision, 2026-09). \"backfill\":",
                    "hedges is offered as a normal candidate predictor -- its pre-2017/2022-2023",
                    "gaps are already filled from the nearest real year upstream in",
                    "dataPrep_Monitor (loadCovariates()/loadHabitatCovariates()/",
                    "occurrencePrepGerHabitat()), so this doesn't invent new fill logic, it",
                    "just re-exposes an already-backfilled column. See covariatePredictorColumns()."),

    ## Spatial blocking parameters --------------------------------------------------
    defineParameter("kFolds", "numeric", 5, NA, NA,
                    "Number of CV folds/blocks."),
    defineParameter("maxBlockSizeEuropeM", "numeric", 1500000, NA, NA,
                    "Maximum spatial block size (m) for the European climate SDM."),
    defineParameter("minBlockSizeEuropeM", "numeric", 200000, NA, NA,
                    "Minimum spatial block size (m) for the European climate SDM."),
    defineParameter("maxBlockSizeGerHabitatM", "numeric", 200000, NA, NA,
                    "Maximum spatial block size (m) for the German habitat SDM."),
    defineParameter("maxBlockSizeGerLandscapeM", "numeric", 200000, NA, NA,
                    "Maximum spatial block size (m) for the German landscape SDM."),
    defineParameter("blockSizeFloorMultiplier", "numeric", 2, NA, NA,
                    "Multiplier applied to a scale's covariate resolution to floor its",
                    "autocorrelation-derived spatial block size (see determineBlockSize()):",
                    "below resolution x this multiplier, adjacent points can share a covariate",
                    "cell across train/test folds. Wiedenroth et al./Lisa Hildebrand's v2 code use",
                    "2x at every scale EXCEPT the 30km landscape variant, which uses 1x -- i.e. even",
                    "the source methodology treats this as scale-dependent, not a fixed constant.",
                    "Exposed as one shared multiplier (not hardcoded per call site) so it can be",
                    "revisited without touching spatialBlockingGerHabitat()/GerLandscape()."),

    ## Collinearity parameters -------------------------------------------------------
    defineParameter("collinearityThreshold", "numeric", 0.7, NA, NA,
                    "Absolute Spearman correlation threshold above which one of a pair of",
                    "predictors is dropped."),
    defineParameter("collinearityUnivar", "character", "gam", NA, NA,
                    "Initial univariate model form for block-CV importance ranking",
                    "(auto-refined internally based on data)."),

    ## Bioclim reference (must match dataPrep_Monitor's values exactly) ----------------
    defineParameter("ebba2TrainingYear", "numeric", 2017, NA, NA,
                    "Target year whose bioclim window was used to train the European EBBA2",
                    "climate SDM in dataPrep_Monitor -- used here to deterministically locate",
                    "the matching bioclim_<start>-<end>.tif as the spatial-blocking reference",
                    "grid, instead of guessing from file modification times. Must match",
                    "dataPrep_Monitor's ebba2TrainingYear."),
    defineParameter("climateWindowLength", "numeric", 6, NA, NA,
                    "Rolling window length (years) used to compute the bioclim climatology.",
                    "Must match dataPrep_Monitor's climateWindowLength."),

    ## Scale resolutions (must match dataPrep_Monitor's/models_Monitor's copies) -------
    defineParameter("climateResolutionM", "numeric", 50000, NA, NA,
                    "Resolution (m) of the climate scale -- used, via scaleLabel(), to",
                    "locate that scale's processed covariates and name its inputs/outputs",
                    "subfolders. Must match dataPrep_Monitor's climateResolutionM."),
    defineParameter("habitatResolutionM", "numeric", 200, NA, NA,
                    "Resolution (m) of the habitat scale -- used, via scaleLabel(), to",
                    "locate that scale's processed covariates and name its inputs/outputs",
                    "subfolders. Must match dataPrep_Monitor's habitatResolutionM."),
    defineParameter("landscapeResolutionM", "numeric", 1000, NA, NA,
                    "Resolution (m) of the landscape scale -- used, via scaleLabel(), to",
                    "locate that scale's processed covariates and name its inputs/outputs",
                    "subfolders. Must match dataPrep_Monitor's landscapeResolutionM."),

    ## Species -------------------------------------------------------------------------
    defineParameter("species", "character",
                    c("Vanellus vanellus", "Milvus milvus", "Lanius collurio",
                      "Lullula arborea", "Alauda arvensis", "Saxicola rubetra",
                      "Emberiza calandra", "Emberiza citrinella", "Buteo buteo",
                      "Sturnus vulgaris", "Perdix perdix"), NA, NA,
                    "Latin names of focal species. Must match dataPrep_Monitor's species param.")
  ),
  inputObjects = bindrows(
    #expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
  ),
  outputObjects = bindrows(
    createsOutput("pooledOccurrence", "list",
                  "List with europe/gerHabitat/gerLandscape named-by-species lists of",
                  "occurrence+covariate data.frames, pooled across years where applicable."),
    createsOutput("spatialBlocks", "list",
                  "List with europe/gerHabitat/gerLandscape named-by-species lists of blocks",
                  "objects (real blockCV::cv_spatial() results, or mimicSpatialBlocks() results",
                  "-- same structure either way)."),
    createsOutput("inputsData", "list",
                  "List with europe/gerHabitat/gerLandscape named-by-species lists, each with",
                  "`data` (the final model-ready table: base ID/coordinate/occurrence/foldID",
                  "columns plus resolved predictor columns) and `predictors` (character vector",
                  "of predictor columns used). This is the input for models_Monitor.")
  )
))

doEvent.inputs_Monitor = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      sim <- scheduleEvent(sim, time(sim), "inputs_Monitor", "spatialBlocking")
      sim <- scheduleEvent(sim, time(sim), "inputs_Monitor", "collinearityCheck")
    },

    spatialBlocking = {
      # ! ----- EDIT BELOW ----- ! #
      occurrenceDir <- file.path(inputPath(sim), "response", "processed")
      europeOcc <- poolOccurrenceEurope(file.path(occurrenceDir, "ornitho"), P(sim)$species)
      habitatOcc <- poolOccurrenceGerHabitat(file.path(occurrenceDir, "MhB"), P(sim)$species)
      landscapeOcc <- poolOccurrenceGerLandscape(file.path(occurrenceDir, "territories"), P(sim)$species)

      predictorsDir <- file.path(inputPath(sim), "predictors", "processed")

      if (isTRUE(P(sim)$runSpatialBlocking)) {
        # Deterministic, not a guess from file mtimes: the exact same
        # bioclim window that occurrencePrepEurope() (dataPrep_Monitor)
        # used to train the EBBA2 climate SDM in the first place.
        windowStart <- P(sim)$ebba2TrainingYear - (P(sim)$climateWindowLength - 1)
        bioclimFile <- file.path(predictorsDir, scaleLabel(P(sim)$climateResolutionM),
                                  paste0("bioclim_", windowStart, "-", P(sim)$ebba2TrainingYear, ".tif"))
        if (!file.exists(bioclimFile)) {
          stop("Expected bioclim training file not found: ", bioclimFile,
               "\nCheck that inputs_Monitor's ebba2TrainingYear/climateWindowLength match ",
               "dataPrep_Monitor's, and that prepareClimateData has run.")
        }

        europeBlocks <- spatialBlockingEurope(
          europeOcc, bioclimFile,
          maxBlockSizeM = P(sim)$maxBlockSizeEuropeM,
          minBlockSizeM = P(sim)$minBlockSizeEuropeM, k = P(sim)$kFolds)
        habitatBlocks <- spatialBlockingGerHabitat(
          habitatOcc, file.path(predictorsDir, scaleLabel(P(sim)$habitatResolutionM),
                                 "solar_radiation_habitat.tif"),
          maxBlockSizeM = P(sim)$maxBlockSizeGerHabitatM,
          minBlockSizeM = P(sim)$blockSizeFloorMultiplier * P(sim)$habitatResolutionM,
          k = P(sim)$kFolds)
        landscapeRefFile <- list.files(file.path(predictorsDir, scaleLabel(P(sim)$landscapeResolutionM)),
                                        pattern = "^landuse_.*\\.tif$", full.names = TRUE)[1]
        landscapeBlocks <- spatialBlockingGerLandscape(
          landscapeOcc, landscapeRefFile,
          maxBlockSizeM = P(sim)$maxBlockSizeGerLandscapeM,
          minBlockSizeM = P(sim)$blockSizeFloorMultiplier * P(sim)$landscapeResolutionM,
          k = P(sim)$kFolds)
      } else {
        message("runSpatialBlocking = FALSE -- using mimicSpatialBlocks() instead of real ",
                "spatial CV blocking.")
        europeBlocks <- lapply(europeOcc, mimicSpatialBlocks, k = P(sim)$kFolds)
        habitatBlocks <- lapply(habitatOcc, mimicSpatialBlocks, k = P(sim)$kFolds)
        landscapeBlocks <- lapply(landscapeOcc, mimicSpatialBlocks, k = P(sim)$kFolds)
      }

      sim$pooledOccurrence <- list(europe = europeOcc, gerHabitat = habitatOcc, gerLandscape = landscapeOcc)
      sim$spatialBlocks <- list(europe = europeBlocks, gerHabitat = habitatBlocks, gerLandscape = landscapeBlocks)
      # ! ----- STOP EDITING ----- ! #
    },

    collinearityCheck = {
      # ! ----- EDIT BELOW ----- ! #
      if (is.null(sim$spatialBlocks)) {
        stop("collinearityCheck requires sim$spatialBlocks -- the spatialBlocking event ",
             "must run first (check your module's scheduleEvent order).")
      }

      # Model-ready per-species tables (data+predictors) are still pipeline
      # INPUT (models_Monitor's input), so they live under inputPath(sim) --
      # never under outputPath(sim), which is reserved for model fitting/
      # prediction results. Corrplots, by contrast, are a diagnostic OUTPUT
      # of this collinearity-selection step, scale-specific, so each goes
      # under that scale's own outputPath(sim) subfolder.
      modelReadyDir <- file.path(inputPath(sim), "model_ready")
      climateLabel <- scaleLabel(P(sim)$climateResolutionM)
      habitatLabel <- scaleLabel(P(sim)$habitatResolutionM)
      landscapeLabel <- scaleLabel(P(sim)$landscapeResolutionM)

      resolveModePerScale <- function(scale) {
        if (is.list(P(sim)$predictorsToUse)) extractScaleExtras(P(sim)$predictorsToUse, scale)
        else P(sim)$predictorsToUse
      }
      resolveTablePerScale <- function(scale) extractScaleExtras(P(sim)$speciesPredictorTable, scale)

      europeResult <- collinearityCheckEurope(
        sim$pooledOccurrence$europe, sim$spatialBlocks$europe,
        predictorsToUse = resolveModePerScale("climate"),
        speciesPredictorTable = resolveTablePerScale("climate"),
        corrplotDir = file.path(outputPath(sim), climateLabel, "corrplots"),
        threshold = P(sim)$collinearityThreshold, univar = P(sim)$collinearityUnivar)
      habitatResult <- collinearityCheckGerHabitat(
        sim$pooledOccurrence$gerHabitat, sim$spatialBlocks$gerHabitat,
        predictorsToUse = resolveModePerScale("habitat"),
        speciesPredictorTable = resolveTablePerScale("habitat"),
        corrplotDir = file.path(outputPath(sim), habitatLabel, "corrplots"),
        threshold = P(sim)$collinearityThreshold, univar = P(sim)$collinearityUnivar,
        hedgesTreatment = P(sim)$hedgesTreatment)
      landscapeResult <- collinearityCheckGerLandscape(
        sim$pooledOccurrence$gerLandscape, sim$spatialBlocks$gerLandscape,
        predictorsToUse = resolveModePerScale("landscape"),
        speciesPredictorTable = resolveTablePerScale("landscape"),
        corrplotDir = file.path(outputPath(sim), landscapeLabel, "corrplots"),
        threshold = P(sim)$collinearityThreshold, univar = P(sim)$collinearityUnivar,
        hedgesTreatment = P(sim)$hedgesTreatment)

      sim$inputsData <- list(europe = europeResult, gerHabitat = habitatResult, gerLandscape = landscapeResult)

      # Persist to disk, same shape regardless of the strategy used, so
      # models_Monitor (or a standalone cluster task) can read from
      # inputPath(sim)/model_ready/<scaleLabel>/ directly.
      scaleDirByName <- c(europe = climateLabel, gerHabitat = habitatLabel, gerLandscape = landscapeLabel)
      for (scaleName in names(sim$inputsData)) {
        scaleDir <- file.path(modelReadyDir, scaleDirByName[[scaleName]])
        dir.create(scaleDir, recursive = TRUE, showWarnings = FALSE)
        for (sp in names(sim$inputsData[[scaleName]])) {
          spClean <- gsub(" ", "_", sp)
          saveRDS(sim$inputsData[[scaleName]][[sp]]$data,
                  file.path(scaleDir, paste0(spClean, "_inputs.rds")))
          saveRDS(sim$inputsData[[scaleName]][[sp]]$predictors,
                  file.path(scaleDir, paste0(spClean, "_predictors.rds")))
        }
      }
      # ! ----- STOP EDITING ----- ! #
    },

    warning(noEventWarning(sim))
  )
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # ! ----- EDIT BELOW ----- ! #

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
