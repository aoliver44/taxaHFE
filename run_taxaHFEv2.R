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
    taxaHFE.R [--subject_identifier=<subject_colname> --label=<label> --feature_type=<feature_type> --sample_fraction=<proportion> --standardized=<TRUE/FALSE> --abundance=<decimal> --prevalence=<decimal> --cor_level=<correlation_level> --format_metaphlan=<TRUE/FALSE> --write_old_files=<TRUE/FALSE> --lowest_level=<integer> --ncores=<ncores>] <input_metadata> <input> <output>
    
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
    --cor_level level of initial correlation filter [default: 0.95]
    --format_metaphlan tells program to expect the desired hData style format, otherwise it attempts to coerce into format [default: FALSE]
    --write_old_files write individual level files and old HFE files [default: TRUE]
    --ncores number of cpu cores to use [default: 2]
Arguments:
    input_meta path to metadata input (txt | tsv | csv)
    input path to input file from hierarchical data (i.e. hData data) (txt | tsv | csv)
    output output file name (csv)

' -> doc

opt <- docopt::docopt(doc, version = 
                        'taxaHFE.R v2\n\n')

## load functions ==============================================================

source("/home/docker/tree.R")

## arg tests ===================================================================
# opt <- data.frame(subject_identifier = character(),
#                   label = character(),
#                   feature_type = character(),
#                   ncores = numeric(),
#                   cor_level = numeric(),
#                   sample_fraction = numeric(),
#                   abundance = numeric(),
#                   prevalence = numeric(),
#                   standardized = character(),
#                   write_old_files = character(),
#                   format_metaphlan=character(),
#                   lowest_level = numeric(),
#                   input_metadata = character(),
#                   input = character(),
#                   output = character())
# opt <- opt %>% tibble::add_row(
#   subject_identifier = "subject_id",
#   label = "cluster",
#   feature_type = "factor",
#   write_old_files = "TRUE",
#   abundance = 0.0001,
#   prevalence = 0.01,
#   standardized = "TRUE",
#   sample_fraction = 1,
#   cor_level = 0.95,
#   format_metaphlan = "TRUE",
#   lowest_level = 2,
#   ncores = 4,
#   input_metadata = "/home/docker/example_inputs/metadata.txt",
#   input = "/home/docker/example_inputs/microbiome_data.txt",
#   output = "/home/docker/example_inputs/output.csv"
# )

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
hData <- read_in_microbiome(input = opt$input, meta = metadata, abundance = opt$abundance, format_metaphlan = opt$format_metaphlan, cores = opt$ncores)

## Build tree ==================================================================
cat("\n\n", "###########################\n", "Building Tree...\n", "###########################\n\n")
cat("This may take a few minutes depending on how many features you have.\n")
hTree <- build_tree(hData, filter_prevalence = opt$prevalence, filter_mean_abundance = opt$abundance)

## Main competition ============================================================

cat("\n\n", "###########################\n", "Competing Tree...\n", "###########################\n\n")

competed_tree <- compete_tree(
  hTree,
  lowest_level = opt$lowest_level,
  max_depth = 1000,
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
flattened_df <- flatten_tree_with_metadata(competed_tree)
colnames(flattened_df) <- gsub(pattern = "abundance\\.", replacement = "", x = colnames(flattened_df))

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

output_nosf <- merge(metadata, output_nosf, by = "subject_id")
output_nosf <- output_nosf %>% janitor::clean_names()
readr::write_delim(file = paste0(tools::file_path_sans_ext(opt$output), "_no_sf.csv"), x = output_nosf, delim = ",")

## Final filter SF =============================================================
cat("\n\n", "###########################\n", "Final Filter...\n", "###########################\n\n")

output_sf <- rf_competition_sf(df = output_nosf, metadata = metadata, 
                  feature_of_interest = "feature_of_interest", 
                  subject_identifier = "subject_id", feature_type = opt$feature_type, 
                  ncores = opt$ncores, nperm = (nperm + 230), 
                  sample_fraction = 
                    calc_class_frequencies(metadata, opt$feature_type, 
                                           feature = "feature_of_interest", 
                                           sample_fraction = opt$sample_fraction), 
                  output = opt$output)

readr::write_delim(file = opt$output, x = output_sf, delim = ",")
cat("\n\n Features (no super filter): ", (ncol(output_nosf) - 2))
cat("\n Features (super filter): ", (NCOL(output_sf) - 2), "\n\n")

## write old files  ============================================================
if (opt$write_old_files == TRUE) {
  cat("\n\n", "###########################\n", "Writing old files...\n", "###########################\n\n")
  
  write_summary_files(input = flattened_df, metadata = metadata, output = opt$output)
  write_old_hfe(input = flattened_df, output = opt$output)
}

save.image(file = paste0(tools::file_path_sans_ext(opt$output), ".RData"), safe = TRUE)
