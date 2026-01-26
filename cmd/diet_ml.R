#!/usr/bin/env Rscript

source("lib/options.R")
source("lib/diet_ml_funcs.R")
source("lib/shap_funcs.R")
source("lib/tree.R")

# load flags
# to use this code line-by-line in the Rstudio context, commandArgs can be overloaded to specify the desired flags
# Example cmd:
# command <- "example_inputs/bike_share_day.csv -o test_outputs -s instant -l cnt --model ridge -t numeric --metric rsq --tune_time 1 --seed 1234 --shap -n 2"
# commandArgs <- function(x) { unlist(strsplit(command, split = " ")) }

# these will be used by the argparser
opts <- load_diet_ml_args()
program <- paste0("dietml")

## initiate logger
initiate_logger(opts_object = opts, program = program)
logger::log_info("Command seen: {paste(commandArgs(), collapse = ' ')}")

## read in data
data <- read_in_metadata(input = opts$DATA,
                             subject_identifier = opts$subject_identifier,
                             label = opts$label, 
                             feature_type = opts$feature_type, 
                             random_effects = FALSE, 
                             limit_covariates = FALSE, 
                             k = NULL,
                             cores = (opts$ncores * opts$parallel_workers))

## split data
tr_te_split <- rsample::initial_split(data, prop = as.numeric(opts$train_split), strata = feature_of_interest)
train_data <- rsample::training(tr_te_split)
test_data  <- rsample::testing(tr_te_split)

## run dietML

## write dietML objects to file (if people want the output files that
readr::write_csv(x = as.data.frame(train_data), 
                 file = paste0(opts$output_dir, "/", program, "_train_", opts$seed, ".csv"), 
                 append = FALSE)
readr::write_csv(x = as.data.frame(test_data), 
                 file = paste0(opts$output_dir, "/", program, "_test_", opts$seed, ".csv"), 
                 append = FALSE)
logger::log_info("Training and testing data written to file")

## pass to dietML if selected
shap_inputs <- run_dietML(train = as.data.frame(train_data), 
                          test = as.data.frame(test_data), 
                          model = opts$model,
                          program = program, 
                          seed = opts$seed, 
                          random_effects = FALSE, 
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
                train = as.data.frame(train_data), 
                test = as.data.frame(test_data), 
                feature_type = opts$feature_type, 
                parallel_workers = opts$parallel_workers
  )
}
