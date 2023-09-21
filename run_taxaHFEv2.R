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

## load libraries =====================================================
# source("/home/docker/tree.R")
source("/scripts/utilities/tree.R")

## add commandline options =====================================================

'Hierarchical feature engineering (HFE) for the reduction of features with respects to a factor or regressor
Usage:
    taxaHFE.R [options] <METADATA> <DATA> <OUTPUT>
    
Options:
    -h --help                         Show help text.
    -v --version                      Show version.
    -s --subject_identifier <string>  Metadata column name containing subject IDs [default: subject_id]
    -l --label <string>               Metadata column name of interest for ML [default: cluster]
    -t --feature_type <string>        Is the ML label a factor or numeric [default: factor]
    -f --sample_fraction <float>      Only let rf see a fraction of total data [default: 1]
    -a --abundance <float>            Minimum mean abundance of feature [default: 0.0001]
    -p --prevalence <float>           Minimum prevalence of feature [default: 0.01]
    -L --lowest_level <int>           Most general level allowed to compete [default: 2]
    -m --max_depth <int>              How many hierarchical levels should be allowed to compete [default: 1000]
    -c --cor_level <float>            Initial pearson correlation filter [default: 0.95]
    -n --ncores <int>                 Number of cpu cores to use [default: 2]
    -d --disable_super_filter         Disable running of the super filter (final forest competition)
    -w --write_old_files              Write individual level files and old HFE files
    -W --write_flattened_tree         Write a compressed backup of the entire competed tree
    -D --write_both_outputs           Write an output for pre and post super filter results, overridden by --disable_super_filter
    --nperm <int>                     Number of RF permutations [default: 40]
    --seed <numeric>                  Set a random numeric seed, default is to use system time

Arguments:
    METADATA path to metadata input (txt | tsv | csv)
    DATA path to input file from hierarchical data (i.e. hData data) (txt | tsv | csv)
    OUTPUT output file name (csv)

' -> doc

# these options will be converted to numeric by load_docopt
numeric_options <- c("sample_fraction", "abundance", "prevalence", "lowest_level", "max_depth", "cor_level", "ncores", "nperm")
# to use this code line-by-line in the Rstudio context, commandArgs can be overloaded to specify the desired flags
# ex. commandArgs <- function(x) { "-s Sample -l Category -L 3 example_inputs/metadata.txt example_inputs/microbiome_data.txt example_inputs/out.txt" }
# these will be used by the options loader
opt <- load_docopt(doc, version = 'taxaHFE.R v2.0\n\n', to_convert = numeric_options)

## Run main ====================================================================

## set random seed
set_seed_func(opt$seed)

## parameters specified
cat("\n","Parameters specified: \n")
cat(paste0("--subject_identifier: ", opt$subject_identifier), "\n")
cat(paste0("--label: ", opt$label), "\n")
cat(paste0("--feature_type: ", opt$feature_type), "\n")
cat(paste0("--sample_fraction: ", opt$sample_fraction), "\n")
cat(paste0("--abundance: ", opt$abundance), "\n")
cat(paste0("--prevalence: ", opt$prevalence), "\n")
cat(paste0("--lowest_level: ", opt$lowest_level), "\n")
cat(paste0("--max_depth: ", opt$max_depth), "\n")
cat(paste0("--cor_level: ", opt$cor_level), "\n")
cat(paste0("--write_old_files: ", opt$write_old_files), "\n")
cat(paste0("--ncores: ", opt$ncores), "\n")
cat(paste0("--nperm: ", opt$nperm), "\n")
cat(paste0("OUTPUT: ", opt$OUTPUT))

## check for inputs and read in read in =======================================================
cat("\n\n", "###########################\n", "Reading in data...\n", "###########################")

## metadata file
metadata <- read_in_metadata(input = opt$METADATA, 
                             subject_identifier = opt$subject_identifier, 
                             label = opt$label)

## hierarchical data file ==========================================================
hData <- read_in_hierarchical_data(input = opt$DATA, 
                            meta = metadata, 
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
  sample_fraction = calc_class_frequencies(
    input = metadata,
    feature_type = opt$feature_type,
    sample_fraction = opt$sample_fraction
  ),
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
