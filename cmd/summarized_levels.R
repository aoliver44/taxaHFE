#!/usr/bin/env Rscript 

## SCRIPT: summarized_levels.R ===============================================
## DATE:   Jun 26, 2024
##
## PURPOSE: To create the summarized levels and assess their performance in
##          the dietML pipeline

## load libraries & functions ==================================================
source("lib/tree.R")
source("lib/options.R")
source("lib/methods.R")

## add commandline options =====================================================

# to use this code line-by-line in the Rstudio context, commandArgs can be overloaded to specify the desired flags
# ex. commandArgs <- function(x) { c("example_inputs/metadata.txt", "example_inputs/microbiome_data.txt", "example_inputs/out.csv", "-s", "Sample", "-l", "Category", "-L", "3", "-n", "4", "--seed", "42", "--train_split", "0.8") }
# these will be used by the argparser
opts <- load_taxa_hfe_ml_args()

## Run main ====================================================================

## set random seed
set_seed_func(opts$seed)

## set target list for dietML input objects
diet_ml_inputs <- list()

## check for inputs and read in read in ========================================
## metadata file
metadata <- read_in_metadata(input = opts$METADATA,
                             subject_identifier = opts$subject_identifier,
                             label = opts$label, 
                             feature_type = opts$feature_type, 
                             random_effects = opts$random_effects, 
                             limit_covariates = TRUE, 
                             k = opts$k_splits)

## hierarchical data file ======================================================
hData <- read_in_hierarchical_data(input = opts$DATA,
                                   metadata = metadata,
                                   cores = opts$ncores)

## set initial test-train split for ML methods =================================
tr_te_split <- rsample::initial_split(metadata, prop = as.numeric(opts$train_split), strata = feature_of_interest)
train_metadata <- rsample::training(tr_te_split)
test_metadata  <- rsample::testing(tr_te_split)

# Run taxaHFE-ML
method_levels(hdata = hData,
              metadata = metadata,
              prevalence = opts$prevalence,
              abundance = opts$prevalence,
              lowest_level = opts$lowest_level,
              max_level = opts$max_level,
              cor_level = 0.1, # low cor_level for speed - makes almost everything a correlation battle (becase the battles dont matter, only the tree)
              ncores = opts$ncores,
              feature_type = opts$feature_type,
              nperm = opts$nperm,
              disable_super_filter = opts$disable_super_filter,
              write_both_outputs = opts$write_both_outputs,
              write_flattened_tree = opts$write_flattened_tree,
              target_list = diet_ml_inputs,
              col_names = colnames(hData)[2:NCOL(hData)],
              output = opts$OUTPUT,
              seed = opts$seed,
              random_effects = random_effects
)

## make sure test train in each item of list
split_train_data(target_list = diet_ml_inputs, attribute_name = "train_test_attr", seed = opts$seed)

## create df for dietML to parse
diet_ml_input_df <- extract_attributes(items_list = diet_ml_inputs)

## write dietML objects to file (if people want the output files that
## went into dietML)
write_list_to_csv(target_list = diet_ml_inputs, directory = dirname(opts$OUTPUT))

## pass to dietML if selected
run_diet_ml(input_df = diet_ml_input_df, n_repeat = 1)
