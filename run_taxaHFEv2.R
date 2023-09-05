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
    -s --subject_identifier metadata column name containing subject IDs [default: subject_id]
    -l --label metadata column name of interest for ML [default: cluster]
    -t --feature_type is the ML label a factor or numeric [default: factor]
    -f --sample_fraction only let rf see a fraction of total data [default: 1]
    -a --abundance feature abundance filter [default: 0.0001]
    -p --prevalence feature prevalence filter [default: 0.01]
    -L --lowest_level most general level allowed to compete [default: 2]
    -m --max_depth how many hierarchical levels should be allowed to compete [default: 1000]
    -c --cor_level initial pearson correlation filter [default: 0.95]
    -w --write_old_files write individual level files and old HFE files [default: TRUE]
    -n --ncores number of cpu cores to use [default: 2]
    --seed set a random numeric seed, default is to use system time
Arguments:
    METADATA path to metadata input (txt | tsv | csv)
    DATA path to input file from hierarchical data (i.e. hData data) (txt | tsv | csv)
    OUTPUT output file name (csv)

' -> doc

opt <- docopt::docopt(doc, version = 
                        'taxaHFE.R v2.0\n\n')

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
#                   seed = character(),
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
#   lowest_level = 3,
#   max_depth = 1000,
#   ncores = 4,
#   seed = 42,
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

cat("\n", "###########################\n", "Competing Tree...\n", "###########################\n\n")

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
cat("\n", "############################################\n", "Flattening tree and writing final output...\n", "############################################\n\n")

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
  cat("\n", "###########################\n", "Writing old files...\n", "###########################\n\n")
  
  write_summary_files(input = flattened_df, metadata = metadata, output = opt$OUTPUT)
  write_old_hfe(input = flattened_df, output = opt$OUTPUT)
}

save.image(file = paste0(tools::file_path_sans_ext(opt$OUTPUT), ".RData"), safe = TRUE)
