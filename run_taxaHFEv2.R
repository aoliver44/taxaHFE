#!/usr/bin/env Rscript 

## SCRIPT: taxaHFE.R ===============================================
## AUTHOR: Andrew Oliver
## DATE:   Aug 14, 2023
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
    taxaHFE.R [--subject_identifier=<subject_colname> --label=<label> --feature_type=<feature_type> --input_covariates=<path> --sample_fraction=<decimal> --standardized=<TRUE/FALSE> --abundance=<decimal> --prevalence=<decimal> --format_metaphlan=<format> --cor_level=<correlation_level> --write_old_files=<TRUE/FALSE> --lowest_level=<integer> --ncores=<ncores>] <input_metadata> <input> <output>
    
Options:
    -h --help  Show this screen.
    -v --version  Show version.
    --subject_identifier name of columns with subject IDs [default: subject_id]
    --label response feature of interest for classification [default: cluster]
    --feature_type of response i.e. numeric or factor [default: factor]
    --input_covariates path to input covariates [default: FALSE]
    --sample_fraction decimal, only let random forests see a fraction of total data [default: 1]
    --standardized is the sum total feature abundance between subjects equal [default: TRUE]
    --abundance pre taxaHFE abundance filter [default: 0.0001]
    --prevalence pre taxaHFE prevalence filter [default: 0.01]
    --lowest_level is the most general level allowed to compete [default: 3]
    --format_metaphlan tells program to expect the desired hData style format, otherwise it attempts to coerce into format [default: FALSE]
    --cor_level level of initial correlation filter [default: 0.95]
    --write_old_files write individual level files and old HFE files [default: TRUE]
    --ncores number of cpu cores to use [default: 2]
Arguments:
    input_meta path to metadata input (txt | tsv | csv)
    input path to input file from hierarchical data (i.e. hData data) (txt | tsv | csv)
    output output file name (txt)

' -> doc

opt <- docopt::docopt(doc, version = 
                        'taxaHFE.R v2\n\n')

## arg tests ===================================================================
opt <- data.frame(subject_identifier = character(),
                  label = character(),
                  feature_type = character(),
                  ncores = numeric(),
                  format_metaphlan = character(),
                  cor_level = numeric(),
                  sample_fraction = numeric(),
                  abundance = numeric(),
                  prevalence = numeric(),
                  standardized = character(),
                  write_old_files = character(),
                  lowest_level = numeric(),
                  input_covariates = character(),
                  input_metadata = character(),
                  input = character(),
                  output = character())
opt <- opt %>% tibble::add_row(
  subject_identifier = "Sample",
  label = "Study.Group",
  feature_type = "factor",
  format_metaphlan = "FALSE",
  write_old_files = "FALSE",
  abundance = 0.0001,
  prevalence = 0.01,
  standardized = "FALSE",
  sample_fraction = 1,
  cor_level = 0.95,
  lowest_level = 2,
  ncores = 4,
  input_metadata = "/home/docker/example_inputs/ERAWIJANTARI_metadata.tsv",
  input_covariates = "FALSE",
  input = "/home/docker/example_inputs/ERAWIJANTARI_species.tsv",
  output = "/home/docker/example_inputs/out.txt"
)

## load functions ==============================================================

source("/home/docker/tree.R")

## Run main ====================================================================

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

## read in metadata file =======================================================
## rename the subject_identifier to subject_id and
## rename the label to feature_of_interest
## metadata, should be in tab or comma separated format

metadata <- read_in_metadata(input = opt$input_metadata, 
                             subject_identifier = opt$subject_identifier, 
                             label = opt$label)

## read in microbiome ==========================================================
## read in data, should be in tab or comma separated format

hData <- read_in_microbiome(input = opt$input, meta = metadata, cores = opt$ncores)

## write old files =============================================================
#write_summary_files(input = hData, output = opt$output)

## Build tree ==================================================================
hTree <- build_tree(hData, filter_prevalence = opt$prevalence, filter_mean_abundance = opt$abundance)

## Main competition ============================================================
competed_tree <- compete_tree(
  hTree,
  lowest_level = opt$lowest_level,
  corr_threshold = opt$cor_level,
  metadata = metadata,
  ncores = opt$ncores,
  feature_type = opt$feature_type,
  nperm = nperm,
  sample_fraction = calc_class_frequencies(
    input = metadata,
    feature_type = opt$feature_type,
    sample_fraction = opt$sample_fraction
  ),
)

## 
# Flatten the tree and tree decisions
flattened_df <- flatten_tree_with_metadata(competed_tree)
# filter to only winners
flattened_df_winners <- flattened_df %>% dplyr::filter(., winner == TRUE)

## write output
output_nosf <- flattened_df_winners %>%
  dplyr::select(., name, 10:dplyr::last_col()) %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames(., var = "name") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "subject_id")

output_nosf$subject_id <- gsub(pattern = "abundance\\.", replacement = "", x = output_nosf$subject_id)
output_nosf <- merge(metadata, output_nosf, by = "subject_id")
readr::write_delim(file = paste0(tools::file_path_sans_ext(opt$output), "taxaHFEv2_no_sf.csv"), x = output_nosf, delim = ",")

## SF 

