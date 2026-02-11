#!/usr/bin/env Rscript

## SCRIPT: taxaHFE-ML.R ===============================================
## DATE:   Jun 26, 2024
##
## PURPOSE: To compress feature space of hierarchical organized data while 
##          considering data leakage. This will assess the performance (in a
##          ML model) of taxaHFE feature reduction

## load libraries & functions ==================================================
source("lib/requirements.R")
source("lib/tree.R")
source("lib/diet_ml_funcs.R")
source("lib/shap_funcs.R")
source("lib/options.R")
source("lib/methods.R")

## add commandline options =====================================================

# to use this code line-by-line in the Rstudio context, commandArgs can be overloaded to specify the desired flags
# ex. command <- "example_inputs/metadata.txt example_inputs/microbiome_data.txt -o test_outputs -s Sample -l Category --seed 1234 --shap -n 2"
# commandArgs <- function(x) { unlist(strsplit(command, split = " ")) }
# these will be used by the argparser

# these will be used by the argparser
opts <- load_taxa_hfe_ml_args()
program <- paste0("taxahfe_ml", ifelse(opts$disable_super_filter, "_no_sf", "_sf"))

## initiate logger
initiate_logger(opts_object = opts, program = program)
logger::log_info("Command seen: {paste(commandArgs(), collapse = ' ')}")

## Run main ====================================================================

## check for inputs and read in read in ========================================
## metadata file
metadata <- read_in_metadata(input = opts$METADATA,
                             subject_identifier = opts$subject_identifier,
                             label = opts$label, 
                             feature_type = opts$feature_type, 
                             random_effects = opts$random_effects, 
                             limit_covariates = TRUE, 
                             k = opts$k_splits,
                             cores = (opts$ncores * opts$parallel_workers))

## hierarchical data file ======================================================
hierarchical_data <- read_in_hierarchical_data(input = opts$DATA,
                                   metadata = metadata,
                                   cores = (opts$ncores * opts$parallel_workers))

## set initial test-train split for ML methods =================================
tr_te_split <- rsample::initial_split(metadata, prop = as.numeric(opts$train_split), strata = feature_of_interest)
train_metadata <- rsample::training(tr_te_split)
test_metadata  <- rsample::testing(tr_te_split)
  
# Run taxaHFE-ML
diet_ml_inputs <- method_taxa_hfe_ml(
  h_data = hierarchical_data,
  metadata = metadata,
  prevalence = opts$prevalence,
  abundance = opts$abundance,
  lowest_level = opts$lowest_level,
  max_level = opts$max_level,
  cor_level = opts$cor_level,
  vif_threshold = opts$vif_threshold,
  ncores = (opts$ncores * opts$parallel_workers),
  feature_type = opts$feature_type,
  nperm = opts$nperm,
  disable_super_filter = opts$disable_super_filter,
  train_metadata = train_metadata,
  test_metadata = test_metadata,
  seed = opts$seed,
  random_effects = opts$random_effects,
  write_flattened_tree = opts$write_flattened_tree,
  output = opts$output_dir
)

## write dietML objects to file (if people want the output files that
readr::write_csv(x = as.data.frame(diet_ml_inputs[1]), 
                 file = paste0(opts$output_dir, "/", program, "_train_", opts$seed, ".csv"), 
                 append = FALSE)
readr::write_csv(x = as.data.frame(diet_ml_inputs[2]), 
                 file = paste0(opts$output_dir, "/", program, "_test_", opts$seed, ".csv"), 
                 append = FALSE)

## pass to dietML if selected
shap_inputs <- run_dietML(train = as.data.frame(diet_ml_inputs[1]), 
                          test = as.data.frame(diet_ml_inputs[2]), 
                          model = opts$model,
                          program = program, 
                          seed = opts$seed, 
                          random_effects = opts$random_effects, 
                          folds = opts$folds, 
                          cv_repeats = opts$cv_repeats, 
                          ncores = opts$ncores, 
                          parallel_workers = opts$parallel_workers, 
                          tune_length = opts$tune_length, 
                          tune_stop = opts$tune_stop, 
                          tune_time = opts$tune_time, 
                          metric = opts$metric, 
                          label = opts$label, 
                          output = opts$output, 
                          feature_type = opts$feature_type, 
                          shap = opts$shap, 
                          cor_level = opts$cor_level, 
                          info_gain_n = opts$info_gain_n
)

## run shap analysis if requested
if (opts$shap) {
  shap_analysis(label = opts$label, 
                output = opts$output, 
                model = opts$model, 
                filename = paste0(program, "_", opts$seed), 
                shap_inputs = shap_inputs, 
                train = as.data.frame(diet_ml_inputs[1]), 
                test = as.data.frame(diet_ml_inputs[2]), 
                feature_type = opts$feature_type, 
                parallel_workers = opts$parallel_workers
  )
}

# also run summarized levels if requested
diet_ml_inputs_levels <- list()
if (opts$summarized_levels) {
  diet_ml_inputs_levels <- method_levels(
    h_data = hierarchical_data, 
    metadata = metadata, 
    prevalence = opts$prevalence,
    abundance = opts$abundance,
    ncores = opts$ncores,
    disable_super_filter = opts$disable_super_filter,
    seed = opts$seed
  )
  
  for (level in seq(1:length(diet_ml_inputs_levels))) {
    ## save them to file
    readr::write_csv(x =diet_ml_inputs_levels[level][[1]], 
                     file = paste0(opts$output_dir, "/", "summarized_level_", level, "_", opts$seed, ".csv"), 
                     append = FALSE)
    
    ## seperate out train and test for summarized levels, based on the 
    ## original splits at the top
    train_levels <- diet_ml_inputs_levels[level][[1]] %>% 
      dplyr::filter(., subject_id %in% train_metadata$subject_id)
    test_levels <- diet_ml_inputs_levels[level][[1]] %>% 
      dplyr::filter(., subject_id %in% test_metadata$subject_id)
    
    ## feed train and test to dietml
    shap_inputs <- run_dietML(train = train_levels, 
                              test = test_levels, 
                              model = opts$model,
                              program = paste0("summarized_level_", level), 
                              seed = opts$seed, 
                              random_effects = opts$random_effects, 
                              folds = opts$folds, 
                              cv_repeats = opts$cv_repeats, 
                              ncores = opts$ncores, 
                              parallel_workers = opts$parallel_workers, 
                              tune_length = opts$tune_length, 
                              tune_stop = opts$tune_stop, 
                              tune_time = opts$tune_time, 
                              metric = opts$metric, 
                              label = opts$label, 
                              output = opts$output, 
                              feature_type = opts$feature_type, 
                              shap = opts$shap, 
                              cor_level = opts$cor_level, 
                              info_gain_n = opts$info_gain_n
    )
    
    if (opts$shap) {
      shap_analysis(label = opts$label, 
                    output = opts$output, 
                    model = opts$model, 
                    filename = paste0("summarized_level_", level, "_", opts$seed), 
                    shap_inputs = shap_inputs, 
                    train = as.data.frame(train_levels), 
                    test = as.data.frame(test_levels), 
                    feature_type = opts$feature_type, 
                    parallel_workers = opts$parallel_workers
      )
    }
  }
}

