#!/usr/bin/env Rscript 

## SCRIPT: taxaHFE.R ===============================================
## AUTHOR: Andrew Oliver
## DATE:   Aug 14, 2023
##
## PURPOSE: To compress feature space of hierarchical organized data

## docker info =================================================================

## docker command:
## to do: are we going to eventually make the docker image have an entrypoint?
#docker run --rm -v `PWD`/:/home/docker -w /home/docker aoliver44/taxa_hfe:latest

## set working dir to /home for the docker container
setwd("/home/docker")

## load libraries & functions ==================================================
source("/scripts/tree.R")
source("/scripts/options.R")

## add commandline options =====================================================

# to use this code line-by-line in the Rstudio context, commandArgs can be overloaded to specify the desired flags
# ex. commandArgs <- function(x) { c("example_inputs/metadata.txt", "example_inputs/microbiome_data.txt", "example_inputs/out.txt", "-s", "Sample", "-l", "Category", "-L", "3", "-n", "4", "--seed", "42") }
# these will be used by the argparser
opt <- load_args('taxaHFE.R v2.11', 1)

## Run main ====================================================================

## set random seed
set_seed_func(opt$seed)

## parameters specified
cat("\n","Parameters specified: \n")
cat(paste0("--subject_identifier: ", opt$subject_identifier), "\n")
cat(paste0("--label: ", opt$label), "\n")
cat(paste0("--feature_type: ", opt$feature_type), "\n")
cat(paste0("--abundance: ", opt$abundance), "\n")
cat(paste0("--prevalence: ", opt$prevalence), "\n")
cat(paste0("--lowest_level: ", opt$lowest_level), "\n")
cat(paste0("--max_depth: ", opt$max_depth), "\n")
cat(paste0("--cor_level: ", opt$cor_level), "\n")
cat(paste0("--write_old_files: ", opt$write_old_files), "\n")
cat(paste0("--ncores: ", opt$ncores), "\n")
cat(paste0("--nperm: ", opt$nperm), "\n")
cat(paste0("--seed: ", opt$seed), "\n")
cat(paste0("OUTPUT: ", opt$OUTPUT))

## check for inputs and read in read in =======================================================
cat("\n\n", "###########################\n", "Reading in data...\n", "###########################")

## metadata file
metadata <- read_in_metadata(input = opt$METADATA, 
                             subject_identifier = opt$subject_identifier, 
                             label = opt$label)

## hierarchical data file ==========================================================
hData <- read_in_hierarchical_data(input = opt$DATA, 
                                   metadata = metadata, 
                                   cores = opt$ncores)

## Build tree ==================================================================
cat("\n\n", "###########################\n", "Building Tree...\n", "###########################\n\n")
cat("This may take a few minutes depending on how many features you have.\n")
hTree <- build_tree(hData, 
                    filter_prevalence = opt$prevalence,
                    filter_mean_abundance = opt$abundance)

## Main competition ============================================================
cat("\n", "###########################\n", "Competing Tree...\n", "###########################\n\n")

competed_tree <- compete_tree(
  hTree,
  lowest_level = opt$lowest_level,
  max_depth = opt$max_depth, # allows for all levels to be competed. Change to 1 for pairwise comparisons
  col_names = colnames(hData)[2:NCOL(hData)],
  corr_threshold = opt$cor_level,
  metadata = metadata,
  ncores = opt$ncores,
  feature_type = opt$feature_type,
  nperm = opt$nperm,
  disable_super_filter = opt$disable_super_filter
)

## Extract information from tree  ==============================================
# Flatten the tree and tree decisions
cat("\n", "############################################\n", "Flattening tree and writing final output...\n", "############################################\n\n")

generate_outputs(
  competed_tree,
  metadata,
  colnames(hData)[2:NCOL(hData)],
  opt$OUTPUT, opt$disable_super_filter,
  opt$write_both_outputs,
  opt$write_old_files,
  opt$write_flattened_tree,
  opt$ncores
)

## message to user finish
cat(" Outputs written! TaxaHFE completed. \n")
