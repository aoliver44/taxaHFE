#!/usr/bin/env Rscript 

## SCRIPT: taxaHFE.R ===============================================
## AUTHOR: Andrew Oliver
## DATE:   Jun 26, 2024
##
## PURPOSE: To compress feature space of hierarchical organized data, and 
##          compete feature engineering methods

## docker info =================================================================

## docker command:
## to do: are we going to eventually make the docker image have an entrypoint?
#docker run --rm -v `PWD`/:/home/docker -w /home/docker aoliver44/taxa_hfe:latest

## set working dir to /home for the docker container
setwd("/home/docker")

## load libraries & functions ==================================================
source("/scripts/tree.R")
source("/scripts/options.R")
source("/scripts/methods.R")

## add commandline options =====================================================

# to use this code line-by-line in the Rstudio context, commandArgs can be overloaded to specify the desired flags
# ex. commandArgs <- function(x) { c("example_inputs/metadata.txt", "example_inputs/microbiome_data.txt", "example_inputs/out.txt", "-s", "Sample", "-l", "Category", "-L", "3", "-n", "4", "--seed", "42") }
# these will be used by the argparser
opt <- load_args("taxaHFE.R v2.11", 2)

## Run main ====================================================================

## set target list for dietML input objects
dietML_inputs <- list()

## set random seed
set_seed_func(opt$seed)

## check for inputs and read in read in ========================================
## metadata file
metadata <- read_in_metadata(input = opt$METADATA,
  subject_identifier = opt$subject_identifier,
  label = opt$label)

## hierarchical data file ======================================================
hData <- read_in_hierarchical_data(input = opt$DATA,
  metadata = metadata,
  cores = opt$ncores)

## set initial test-train split for ML methods =================================
tr_te_split <- rsample::initial_split(metadata, prop = as.numeric(opt$train_split), strata = feature_of_interest)
train_metadata <- rsample::training(tr_te_split)
test_metadata  <- rsample::testing(tr_te_split)


## Run TaxaHFE
method_taxaHFE(hdata = hData,
  metadata = metadata,
  prevalence = opt$prevalence,
  abundance = opt$abundance,
  lowest_level = opt$lowest_level,
  max_level = opt$max_level,
  cor_level = opt$cor_level,
  ncores = opt$ncores,
  feature_type = opt$feature_type,
  nperm = opt$nperm,
  disable_super_filter = opt$disable_super_filter,
  write_both_outputs = opt$write_both_outputs,
  write_flattened_tree = opt$write_flattened_tree,
  target_list = dietML_inputs,
  col_names = colnames(hData)[2:NCOL(hData)],
  output = opt$OUTPUT,
  seed = opt$seed
  )

## Run TaxaHFE-ML
method_taxaHFE_ml(hdata = hData,
  metadata = metadata,
  prevalence = opt$prevalence,
  abundance = opt$abundance,
  lowest_level = opt$lowest_level,
  max_level = opt$max_level,
  cor_level = opt$cor_level,
  ncores = opt$ncores,
  feature_type = opt$feature_type,
  nperm = opt$nperm,
  disable_super_filter = opt$disable_super_filter,
  write_both_outputs = opt$write_both_outputs,
  write_flattened_tree = opt$write_flattened_tree,
  train_split = opt$train_split,
  model = opt$model,
  folds = opt$folds,
  metric = opt$metric,
  tune_length = opt$tune_length,
  tune_time = opt$tune_time,
  tune_stop = opt$tune_stop,
  shap = opt$shap,
  target_list = dietML_inputs,
  output = opt$OUTPUT,
  seed = opt$seed
  )

method_levels(hdata = hData,
               metadata = metadata,
               prevalence = opt$prevalence,
               abundance = opt$abundance,
               lowest_level = opt$lowest_level,
               max_level = opt$max_level,
               cor_level = opt$cor_level,
               ncores = opt$ncores,
               feature_type = opt$feature_type,
               nperm = opt$nperm,
               disable_super_filter = opt$disable_super_filter,
               write_both_outputs = opt$write_both_outputs,
               write_flattened_tree = opt$write_flattened_tree,
               target_list = dietML_inputs,
               col_names = colnames(hData)[2:NCOL(hData)],
               output = opt$OUTPUT,
               seed = opt$seed
)


## make sure test train in each item of list

## pass to dietML
