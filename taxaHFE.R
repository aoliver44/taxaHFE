#!/usr/bin/env Rscript 

## SCRIPT: taxaHFE.R ===============================================
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
source("/home/docker/tree.R")
source("/home/docker/options_new.R")
source("/home/docker/methods.R")

## add commandline options =====================================================

# to use this code line-by-line in the Rstudio context, commandArgs can be overloaded to specify the desired flags
# ex. commandArgs <- function(x) { c("example_inputs/metadata.txt", "example_inputs/microbiome_data.txt", "example_inputs/out.txt", "-s", "Sample", "-l", "Category", "-L", "3", "-n", "4", "--seed", "42", "-wWD") }
# these will be used by the argparser
opt <- load_args("taxaHFE", version = "2.3")

## Run main ====================================================================

## set random seed
set_seed_func(opt$seed)

## set target list for dietML input objects
dietML_inputs <- list()

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
