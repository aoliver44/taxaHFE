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

## add commandline options =====================================================

library(docopt)
'Hierarchical feature engineering (HFE) for the reduction of features with respects to a factor or regressor
Usage:
    taxaHFE.R [options] <METADATA> <DATA> <OUTPUT>
    
Options:
    -h --help  Show this screen.
    -v --version  Show version.
    --subject_identifier name of columns with subject IDs [default: subject_id]
    --label response feature of interest for classification [default: cluster]
    --feature_type of response i.e. numeric or factor [default: factor]
    --sample_fraction only let rf see a fraction of total data [default: 1]
    --standardized sum feature abundance across subjects equal [default: TRUE]
    --abundance pre taxaHFE abundance filter [default: 0.0001]
    --prevalence pre taxaHFE prevalence filter [default: 0.01]
    --lowest_level is the most general level allowed to compete [default: 2]
    --max_depth how deep should comparisons be allowed to go [default: 1000]
    --cor_level level of initial correlation filter [default: 0.95]
    --write_old_files write individual level files and old HFE files [default: TRUE]
    --ncores number of cpu cores to use [default: 2]
Arguments:
    METADATA path to metadata input (txt | tsv | csv)
    DATA path to input file from hierarchical data (i.e. hData data) (txt | tsv | csv)
    OUTPUT output file name (csv)

' -> doc

opt <- docopt::docopt(doc, version = 
                        'taxaHFE.R v2\n\n')

## load functions ==============================================================

source("/scripts/utilities/tree.R")

## arg tests ===================================================================
# opt <- data.frame(subject_identifier = character(),
#                   label = character(),
#                   feature_type = character(),
#                   ncores = numeric(),
#                   cor_level = numeric(),
#                   sample_fraction = numeric(),
#                   abundance = numeric(),
#                   prevalence = numeric(),
#                   write_old_files = character(),
#                   lowest_level = numeric(),
#                   max_depth = numeric(),
#                   METADATA = character(),
#                   DATA = character(),
#                   OUTPUT = character())
# opt <- opt %>% tibble::add_row(
#   subject_identifier = "Sample",
#   label = "Category",
#   feature_type = "factor",
#   write_old_files = "TRUE",
#   abundance = 0.0001,
#   prevalence = 0.01,
#   sample_fraction = 1,
#   cor_level = 0.95,
#   lowest_level = 2,
#   max_depth = 1000,
#   ncores = 4,
#   METADATA = "/home/docker/example_inputs/metadata.txt",
#   DATA = "/home/docker/example_inputs/microbiome_data.txt",
#   OUTPUT = "/home/docker/example_inputs/test/output.csv"
# )

## Run main ====================================================================

## check for inputs ============================================================
cat("\n\n", "###########################\n", "Reading in data...\n", "###########################")

## check and see if clean_files directory exists
cat("\n\n","Checking for METADATA")
if (file.exists(opt$METADATA)) {
  cat("\n",paste0("Using ", opt$METADATA, " as input")) 
} else { stop("Metadata input not found.") }

## check for input file (hierarchical data)
cat("\n","Checking for for input...")
if (file.exists(opt$DATA)) {
  cat("\n",paste0("Using ", opt$DATA, " as input")) 
} else { stop("Input not found.") }

## read in metadata file =======================================================

## rename the subject_identifier to subject_id and
## rename the label to feature_of_interest
## metadata, should be in tab or comma separated format
metadata <- read_in_metadata(input = opt$METADATA, 
                             subject_identifier = opt$subject_identifier, 
                             label = opt$label)

## read in microbiome ==========================================================

## read in data, should be in tab or comma separated format
hData <- read_in_microbiome(input = opt$DATA, meta = metadata, abundance = opt$abundance, cores = opt$ncores)

## Build tree ==================================================================
cat("\n\n", "###########################\n", "Building Tree...\n", "###########################\n\n")
cat("This may take a few minutes depending on how many features you have.\n")
hTree <- build_tree(hData, filter_prevalence = opt$prevalence, filter_mean_abundance = opt$abundance)

## Main competition ============================================================

cat("\n\n", "###########################\n", "Competing Tree...\n", "###########################\n\n")

competed_tree <- compete_tree(
  hTree,
  lowest_level = opt$lowest_level,
  max_depth = as.numeric(opt$max_depth), # allows for all levels to be competed. Change to 1 for pairwise comparisons
  col_names = colnames(hData)[2:NCOL(hData)],
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

## Extract information from tree  ==============================================
# Flatten the tree and tree decisions
cat("\n\n", "############################################\n", "Flattening tree and writing final output...\n", "############################################\n\n")

col_names = colnames(hData)[2:NCOL(hData)]
flattened_df <- flatten_tree_with_metadata(node = competed_tree)
colnames(flattened_df)[11:NCOL(flattened_df)] <- col_names

## filter to only winners
flattened_noSF_winners <- flattened_df %>% 
  dplyr::filter(., winner == TRUE)

## clean names in case of duplicates
flattened_noSF_winners$name <- janitor::make_clean_names(flattened_noSF_winners$name)

## write noSF output
output_nosf <- flattened_noSF_winners %>%
  dplyr::select(., name, 11:dplyr::last_col()) %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames(., var = "name") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "subject_id")

output_nosf <- merge(metadata, output_nosf, by = "subject_id")
readr::write_delim(file = paste0(tools::file_path_sans_ext(opt$OUTPUT), "_no_sf.csv"), x = output_nosf, delim = ",")

## write SF output
output_sf <- flattened_noSF_winners %>%
  dplyr::filter(., sf_winner == TRUE) %>%
  dplyr::select(., name, 11:dplyr::last_col()) %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames(., var = "name") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "subject_id")

output_sf <- merge(metadata, output_sf, by = "subject_id")
readr::write_delim(file = opt$OUTPUT, x = output_sf, delim = ",")

## message to user outputs
cat("\n\n", " Outputs written! \n")

cat("\n\n Features (no super filter): ", (ncol(output_nosf) - 2))
cat("\n Features (super filter): ", (NCOL(output_sf) - 2), "\n\n")

## write old files  ============================================================
if (opt$write_old_files == TRUE) {
  cat("\n\n", "###########################\n", "Writing old files...\n", "###########################\n\n")
  
  write_summary_files(input = flattened_df, metadata = metadata, output = opt$OUTPUT)
  write_old_hfe(input = flattened_df, output = opt$OUTPUT)
}

save.image(file = paste0(tools::file_path_sans_ext(opt$OUTPUT), ".RData"), safe = TRUE)
