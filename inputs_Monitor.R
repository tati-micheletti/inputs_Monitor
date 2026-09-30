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
                   "PredictiveEcology/reproducible@development",
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
    defineParameter("speciesPredictorTable", "list", NULL, NA, NA,
                    "The ONLY source of a species' candidate predictors -- no mode selector, ",
                    "every species must be listed here (a species missing here is a hard ",
                    "error at collinearityCheck time, not a silent fallback). ",
                    "Named list: species -> scale -> character vector of predictor names, e.g.: ",
                    "list(\"Buteo buteo\" = list(habitat = c(\"hedges\", \"grassland\"))). ",
                    "(Populated from speciesConfig_predictors.csv when run via runMe.R.)"),
    defineParameter("dropCollinearPredictors", "logical", FALSE, NA, NA,
                    "FALSE (default): use each species' speciesPredictorTable list exactly as ",
                    "given. TRUE: prune that same list for collinearity via real block-CV ",
                    "selection (select07Blockcv()), capped at 1 predictor per 10 occurrences -- ",
                    "narrows what the table says to consider, never substitutes a different ",
                    "candidate set. A module-level technical/algorithmic toggle, not per-species ",
                    "(matches collinearityThreshold/collinearityUnivar's existing pattern)."),
    defineParameter("spatialTermConfig", "list", NULL, NA, NA,
                    "Which species+scale get projected x/y coordinates added as an extra ",
                    "BRT predictor (a spatial trend-surface term meant to absorb residual ",
                    "regional structure -- NOT a formal random effect, see DECISIONS.md). ",
                    "Named list: species -> scale -> TRUE/FALSE. A species+scale left out, ",
                    "or set FALSE, never gets x/y added. Deliberately opt-in per species+scale ",
                    "rather than blanket-applied -- it can just as easily hurt a model ",
                    "(overfitting to historical geography, reduced transportability to future ",
                    "climate-scale predictions) as help it. ",
                    "(Populated from speciesConfig_general.csv's spatial_term column when run ",
                    "via runMe.R, via extractSpatialTermSpecies() -- see that file.)"),
    defineParameter("resolutionConfig", "list", NULL, NA, NA,
                    "Named list: species -> scale -> resolution (m), or NA (falls back to that ",
                    "scale's shared *ResolutionM default). Used to group species by their ",
                    "resolved resolution for spatial blocking -- each distinct group gets its ",
                    "own reference grid and its own block-size floor. (Populated from ",
                    "speciesConfig_general.csv's resolution_m column via ",
                    "extractResolutionConfig() -- see sharedSpeciesConfig.R.)"),
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

    ## Fitting years (must match dataPrep_Monitor's/models_Monitor's copies) -----------
    defineParameter("habitatYears", "numeric", NULL, NA, NA,
                    "NULL (default): pool every <species>_habitat_<year>.rds file found,",
                    "regardless of year (today's original behavior). Otherwise restrict",
                    "poolOccurrenceGerHabitat() to exactly these years -- confirmed necessary",
                    "2026-09-30: without it, stale files left over from an earlier run with a",
                    "different year range get silently pooled alongside the current run's",
                    "fresh files. Should match dataPrep_Monitor's habitatYears."),
    defineParameter("landscapeYears", "numeric", NULL, NA, NA,
                    "Same rationale as habitatYears, for poolOccurrenceGerLandscape()'s",
                    "<species>_landscape_<year>.rds files. Should match dataPrep_Monitor's",
                    "landscapeYears."),

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
    defineParameter("species", "character", NA_character_, NA, NA,
                    "Latin names of focal species -- no default (errors if unset); supply",
                    "sharedSpecies from sharedConfig.R (repo root), same as",
                    "dataPrep_Monitor's species param, so the two can never silently",
                    "drift apart.")
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
      if (identical(P(sim)$species, NA_character_)) {
        stop("inputs_Monitor's species parameter must be supplied explicitly ",
             "(e.g. sharedSpecies from sharedConfig.R) -- no default roster.")
      }
      sim <- scheduleEvent(sim, time(sim), "inputs_Monitor", "spatialBlocking")
      sim <- scheduleEvent(sim, time(sim), "inputs_Monitor", "collinearityCheck")
    },

    spatialBlocking = {
      # ! ----- EDIT BELOW ----- ! #
      occurrenceDir <- file.path(inputPath(sim), "response", "processed")
      europeOcc <- poolOccurrenceEurope(file.path(occurrenceDir, "ornitho"), P(sim)$species)
      habitatOcc <- poolOccurrenceGerHabitat(file.path(occurrenceDir, "MhB"), P(sim)$species,
                                              years = P(sim)$habitatYears)
      landscapeOcc <- poolOccurrenceGerLandscape(file.path(occurrenceDir, "territories"), P(sim)$species,
                                                   years = P(sim)$landscapeYears)

      predictorsDir <- file.path(inputPath(sim), "predictors", "processed")

      if (isTRUE(P(sim)$runSpatialBlocking)) {
        # Group species by their own resolved resolution per scale (see
        # DECISIONS.md, 2026-09-28) -- each distinct group gets its own
        # reference grid + block-size floor, computed once per group
        # rather than once globally. A species absent from
        # resolutionConfig, or with a blank/NA entry there, falls back to
        # that scale's shared default -- so with nobody overridden (as of
        # today) this produces exactly one group per scale, identical to
        # the previous behaviour.
        groupSpeciesByResolution <- function(species, scale, sharedDefault) {
          resVals <- vapply(species, function(sp) {
            v <- if (!is.null(P(sim)$resolutionConfig)) P(sim)$resolutionConfig[[sp]][[scale]] else NULL
            if (is.null(v) || is.na(v)) sharedDefault else v
          }, numeric(1))
          split(species, resVals)
        }

        # Deterministic, not a guess from file mtimes: the exact same
        # bioclim window that occurrencePrepEurope() (dataPrep_Monitor)
        # used to train the EBBA2 climate SDM in the first place.
        windowStart <- P(sim)$ebba2TrainingYear - (P(sim)$climateWindowLength - 1)

        europeGroups <- groupSpeciesByResolution(names(europeOcc), "climate", P(sim)$climateResolutionM)
        europeBlocks <- list()
        for (resStr in names(europeGroups)) {
          resM <- as.numeric(resStr)
          resLabel <- scaleLabel(resM)
          bioclimFile <- file.path(predictorsDir, resLabel,
                                    paste0("bioclim_", windowStart, "-", P(sim)$ebba2TrainingYear,
                                           "_", resLabel, ".tif"))
          if (!file.exists(bioclimFile)) {
            stop("Expected bioclim training file not found: ", bioclimFile,
                 "\nCheck that inputs_Monitor's ebba2TrainingYear/climateWindowLength match ",
                 "dataPrep_Monitor's, and that prepareClimateData has run.")
          }
          europeBlocks <- c(europeBlocks, spatialBlockingEurope(
            europeOcc[europeGroups[[resStr]]], bioclimFile,
            maxBlockSizeM = P(sim)$maxBlockSizeEuropeM,
            minBlockSizeM = P(sim)$minBlockSizeEuropeM, k = P(sim)$kFolds))
        }

        habitatGroups <- groupSpeciesByResolution(names(habitatOcc), "habitat", P(sim)$habitatResolutionM)
        habitatBlocks <- list()
        for (resStr in names(habitatGroups)) {
          resM <- as.numeric(resStr)
          resLabel <- scaleLabel(resM)
          refPath <- file.path(predictorsDir, resLabel, paste0("solar_radiation_habitat_", resLabel, ".tif"))
          habitatBlocks <- c(habitatBlocks, spatialBlockingGerHabitat(
            habitatOcc[habitatGroups[[resStr]]], refPath,
            maxBlockSizeM = P(sim)$maxBlockSizeGerHabitatM,
            minBlockSizeM = P(sim)$blockSizeFloorMultiplier * resM,
            k = P(sim)$kFolds))
        }

        landscapeGroups <- groupSpeciesByResolution(names(landscapeOcc), "landscape", P(sim)$landscapeResolutionM)
        landscapeBlocks <- list()
        for (resStr in names(landscapeGroups)) {
          resM <- as.numeric(resStr)
          resLabel <- scaleLabel(resM)
          landscapeRefFile <- list.files(file.path(predictorsDir, resLabel),
                                          pattern = "^landuse_.*\\.tif$", full.names = TRUE)[1]
          landscapeBlocks <- c(landscapeBlocks, spatialBlockingGerLandscape(
            landscapeOcc[landscapeGroups[[resStr]]], landscapeRefFile,
            maxBlockSizeM = P(sim)$maxBlockSizeGerLandscapeM,
            minBlockSizeM = P(sim)$blockSizeFloorMultiplier * resM,
            k = P(sim)$kFolds))
        }
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

      resolveTablePerScale <- function(scale) extractScaleExtras(P(sim)$speciesPredictorTable, scale)
      resolveSpatialTermPerScale <- function(scale) extractScaleExtras(P(sim)$spatialTermConfig, scale)

      europeResult <- collinearityCheckEurope(
        sim$pooledOccurrence$europe, sim$spatialBlocks$europe,
        speciesPredictorTable = resolveTablePerScale("climate"),
        dropCollinearPredictors = P(sim)$dropCollinearPredictors,
        spatialTermSpecies = resolveSpatialTermPerScale("climate"),
        corrplotDir = file.path(outputPath(sim), climateLabel, "corrplots"),
        threshold = P(sim)$collinearityThreshold, univar = P(sim)$collinearityUnivar,
        cachePath = cachePath(sim))
      habitatResult <- collinearityCheckGerHabitat(
        sim$pooledOccurrence$gerHabitat, sim$spatialBlocks$gerHabitat,
        speciesPredictorTable = resolveTablePerScale("habitat"),
        dropCollinearPredictors = P(sim)$dropCollinearPredictors,
        spatialTermSpecies = resolveSpatialTermPerScale("habitat"),
        corrplotDir = file.path(outputPath(sim), habitatLabel, "corrplots"),
        threshold = P(sim)$collinearityThreshold, univar = P(sim)$collinearityUnivar,
        cachePath = cachePath(sim))
      landscapeResult <- collinearityCheckGerLandscape(
        sim$pooledOccurrence$gerLandscape, sim$spatialBlocks$gerLandscape,
        speciesPredictorTable = resolveTablePerScale("landscape"),
        dropCollinearPredictors = P(sim)$dropCollinearPredictors,
        spatialTermSpecies = resolveSpatialTermPerScale("landscape"),
        corrplotDir = file.path(outputPath(sim), landscapeLabel, "corrplots"),
        threshold = P(sim)$collinearityThreshold, univar = P(sim)$collinearityUnivar,
        cachePath = cachePath(sim))

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
