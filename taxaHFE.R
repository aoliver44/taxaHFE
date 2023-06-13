#!/usr/bin/env Rscript 

## SCRIPT: taxaHFE.R ===============================================
## AUTHOR: Andrew Oliver
## DATE:   Feb 16, 2022
##
## PURPOSE: To compress feature space of hierarchical organized data

## docker info =================================================================

## docker command:
#docker run --rm -v `PWD`/:/home/docker -w /home/docker aoliver44/taxa_hfe:latest

## set working dir to /home for the docker container
setwd("/home/docker")

## add commandline options =====================================================

library(docopt)
'Hierarchical feature engineering (HFE) for the reduction of features with respects to a factor or regressor
Usage:
    taxaHFE.R [--subject_identifier=<subject_colname> --label=<label> --feature_type=<feature_type> --input_covariates=<path> --subsample=<decimal> --standardized=<TRUE/FALSE> --abundance=<decimal> --prevalence=<decimal> --format_metaphlan=<format> --cor_level=<correlation_level> --write_old_files=<TRUE/FALSE> --ncores=<ncores>] <input_metadata> <input> <output>
    
Options:
    -h --help  Show this screen.
    -v --version  Show version.
    --subject_identifier name of columns with subject IDs [default: subject_id]
    --label response feature of interest for classification [default: cluster]
    --feature_type of response i.e. numeric or factor [default: factor]
    --input_covariates path to input covariates [default: FALSE]
    --subsample decimal, only let random forests see a fraction of total data [default: 1]
    --standardized is the sum total feature abundance between subjects equal [default: TRUE]
    --abundance pre taxaHFE abundance filter [default: 0.0001]
    --prevalence pre taxaHFE prevalence filter [default: 0.01]
    --format_metaphlan tells program to expect the desired hData style format, otherwise it attempts to coerce into format [default: FALSE]
    --cor_level level of initial correlation filter [default: 0.95]
    --write_old_files write individual level files and old HFE files [default: TRUE]
    --ncores number of cpu cores to use [default: 2]
Arguments:
    input_meta path to metadata input (txt | tsv | csv)
    input path to input file from hierarchical data (i.e. hData data) (txt | tsv | csv)
    output output file name (txt)

' -> doc

opt <- docopt::docopt(doc, version = 'taxaHFE.R v1.10\n\n')
#print(opt)
## load libraries ==============================================================

library(dplyr, quietly = T, verbose = F, warn.conflicts = F)
library(janitor, quietly = T, verbose = F, warn.conflicts = F)
library(tidyr, quietly = T, verbose = F, warn.conflicts = F)
library(tibble, quietly = T, verbose = F, warn.conflicts = F)
library(caret, quietly = T, verbose = F, warn.conflicts = F)
library(readr, quietly = T, verbose = F, warn.conflicts = F)
library(reshape2, quietly = T, verbose = F, warn.conflicts = F)
library(ggplot2, quietly = T, verbose = F, warn.conflicts = F)
library(ggsci, quietly = T, verbose = F, warn.conflicts = F)
suppressPackageStartupMessages(library(vegan, quietly = T, verbose = F, warn.conflicts = F))
library(progress, quietly = T, verbose = F, warn.conflicts = F)
library(stringr, quietly = T, verbose = F, warn.conflicts = F)

## set random seed if needed
set.seed(42)
nperm = 10
## helper functions ============================================================

## Negate function ("not in"):
`%!in%` <- Negate(`%in%`)

## suppress warnings
options(warn=-1)

source("/scripts/utilities/taxaHFE_functions.R")

## arg tests ===================================================================
# opt <- data.frame(subject_identifier=character(),
#                   label=character(),
#                   feature_type=character(),
#                   ncores=numeric(),
#                   format_metaphlan=character(),
#                   cor_level=numeric(),
#                   subsample=numeric(),
#                   abundance=numeric(),
#                   prevelance=numeric(),
#                   standardized=character(),
#                   write_old_files=character(),
#                   input_covariates=character(),
#                   input_metadata=character(),
#                   input=character(),
#                   output=character())
# opt <- opt %>% tibble::add_row(
#   subject_identifier = "Sample",
#   label= "Category",
#   feature_type = "factor",
#   format_metaphlan = "TRUE",
#   write_old_files = "FALSE",
#   abundance = 0.0001,
#   prevelance = 0.01,
#   standardized = "TRUE",
#   subsample = 1,
#   cor_level = 0.95,
#   ncores = 4,
#   input_metadata = "/home/docker/example_inputs/metadata.txt",
#   input_covariates = "FALSE",
#   input= "/home/docker/example_inputs/microbiome_data.txt",
#   output = "/home/docker/example_inputs/output.txt"
# )

## check for inputs ============================================================
cat("\n\n", "###########################\n", "Reading in data...\n", "###########################")

## check and see if clean_files directory exists
cat("\n\n","Checking for for input_metadata...")
if (file.exists(opt$input_metadata)) {
  cat("\n",paste0("Using ", opt$input_metadata, " as input")) 
} else { stop("Metadata input not found.") }

## check for input file (hierarchical data)
cat("\n","Checking for for input...")
if (file.exists(opt$input)) {
  cat("\n",paste0("Using ", opt$input, " as input")) 
} else { stop("Input not found.") }

## read in microbiome ==========================================================
## read in data, should be in tab or comma separated format

hData <- read_in_microbiome(input = opt$input)
original_taxa_count <- NROW(hData)

## read in metadata file =======================================================
## rename the subject_identifier to subject_id and
## rename the label to feature_of_interest
## metadata, should be in tab or comma separated format

metadata <- read_in_metadata(input = opt$input_metadata, 
                             subject_identifier = opt$subject_identifier, 
                             label = opt$label)

## run safety checks ===========================================================

mid_safety_checks()

## calculate vector of class frequencies, take 70% of them to leave
## some data out
calc_class_frequencies()

## read in covariates file =====================================================

if (opt$input_covariates != "FALSE") {
  covariates = read_in_covariates(input = opt$input_covariates, subject_identifier = opt$subject_identifier)
}

## if not "metaphlan" format, attempt to convert ===============================

if (opt$format_metaphlan == "FALSE") {
  
  convert_to_hData(input = hData)
  hData <- do.call(rbind, lapply(ls(pattern = "hData_L"), get))
  
}

## Remove very low prevalent features ==========================================

apply_filters(input = hData, 
              standardized = opt$standardized, 
              filter_abundance = opt$abundance, 
              filter_prevalence = opt$prevalence)

## write files for old HFE program =============================================

if (opt$write_old_files == "TRUE") {
  write_old_hfe(input = hData, output = opt$output)
}

## make the dataframe of features that will compete in parent-child competitions.
## Basically this is just the taxonomic (hierarchical) information of all
## the features left at this step.

make_taxa_split_df(input = hData)

## write summarized files and clean up =========================================
## these are the species only, genus only...etc files
## to check against.

if (opt$write_old_files == "TRUE") {
  write_summary_files(input = hData, output = opt$output)
}

## compete! ====================================================================

if (opt$input_covariates == "FALSE") {
  taxaHFE_competition(input = hData, feature_type = opt$feature_type, cores = opt$ncores, output = opt$output)
} else {
  taxaHFE_competition_covariates(input = hData, covariates = covariates, feature_type = opt$feature_type, cores = opt$ncores, output = opt$output)
}
## super filter ================================================================

if (opt$input_covariates == "FALSE") {
  super_filter(input = hData, feature_type = opt$feature_type, cores = opt$ncores, subsample = opt$subsample, output = opt$output)
} else {
  super_filter_covariates(input = hData, covariates = covariates, feature_type = opt$feature_type, cores = opt$ncores, subsample = opt$subsample, output = opt$output)
}

## Write Figure ================================================================

write_figure(input = hData, output = opt$output)


## Write outputs ===============================================================

cat("\n",paste0("Reduced/compressed taxa set from ", original_taxa_count, " taxa to ", (NROW(taxa_only_split_sf)), " (", NROW(taxa_only_split), " if no super filter)\n"))

cat("\n","Output written.  \n\n")
