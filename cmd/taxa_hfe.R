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
# ex. commandArgs <- function(x) { c("metadata_time.txt", "microbiome_time.txt", "out.txt", "-s", "subject_id", "-l", "Intervention", "-L", "3", "-n", "4", "--seed", "42", "-wWD") }
# these will be used by the argparser
opt <- load_args("taxa_hfe")

## Run main ====================================================================

## set random seed
set_seed_func(opt$seed)

## set target list for dietML input objects
diet_ml_inputs <- list()

## check for inputs and read in read in ========================================
## metadata file
metadata <- read_in_metadata(input = opt$METADATA,
                             subject_identifier = opt$subject_identifier,
                             label = opt$label)

## hierarchical data file ======================================================
hData <- read_in_hierarchical_data(input = opt$DATA,
                                   metadata = metadata,
                                   cores = opt$ncores)


# Run TaxaHFE
method_taxa_hfe(hdata = hData,
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
  target_list = diet_ml_inputs,
  col_names = colnames(hData)[2:NCOL(hData)],
  output = opt$OUTPUT,
  seed = opt$seed,
  random_effects = TRUE
  )


## DRAFT CODE - DELETE BEFORE MERGING
hTree <- build_tree(hData,
                    filter_prevalence = 0.01,
                    filter_mean_abundance = 0
)

competed_tree <- compete_tree(
  hTree,
  lowest_level = 3,
  max_level = 1000, # allows for all levels to be competed. Change to 1 for pairwise comparisons
  col_names = colnames(hData)[2:NCOL(hData)],
  corr_threshold = 0.95,
  metadata = metadata,
  ncores = 4,
  feature_type = "factor",
  nperm = 10,
  disable_super_filter = TRUE,
  random_effects = TRUE
)

