#!/usr/bin/env Rscript 

## SCRIPT: taxaHFE.R ===============================================
## DATE:   Jun 26, 2024
##
## PURPOSE: To compress feature space of hierarchical organized data, and 
##          compete feature engineering methods

## load libraries & functions ==================================================
source("lib/tree.R")
source("lib/options.R")
source("lib/methods.R")

## add commandline options =====================================================

# to use this code line-by-line in the Rstudio context, commandArgs can be overloaded to specify the desired flags
# ex. commandArgs <- function(x) { c("example_inputs/metadata.txt", "example_inputs/microbiome_data.txt", "-o", "example_outputs", "-s", "Sample", "-l", "Category", "-L", "3", "-n", "4", "--seed", "42", "-wWD") }
# these will be used by the argparser
opts <- load_taxa_hfe_args()

## Run main ====================================================================

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


# Run TaxaHFE
diet_ml_inputs <- method_taxa_hfe(
  hdata = hData,
  metadata = metadata,
  prevalence = opts$prevalence,
  abundance = opts$abundance,
  lowest_level = opts$lowest_level,
  max_level = opts$max_level,
  cor_level = opts$cor_level,
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
  random_effects = opts$random_effects
)

