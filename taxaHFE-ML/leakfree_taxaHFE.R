#!/usr/bin/env Rscript 

## SCRIPT: leakfree_taxaHFE.R ===============================================
## AUTHOR: Andrew Oliver
## DATE:   April 10, 2024
##
## PURPOSE: To perform data leakage free TaxaHFE and SHAP analysis

## docker info =================================================================

## docker command:
## to do: are we going to eventually make the docker image have an entrypoint?
#docker run --rm -v `pwd`:/home/docker -w /home/docker aoliver44/leakage_free_taxaHFE:latest

## set working dir to /home for the docker container
setwd("/home/docker")

## load libraries & functions ==================================================
source("/scripts/tree.R")
source("/scripts/options.R")

## add commandline options =====================================================

# to use this code line-by-line in the Rstudio context, commandArgs can be overloaded to specify the desired flags
# ex. commandArgs <- function(x) { c("/home/docker/example_inputs/metadata.txt", "/home/docker/example_inputs/microbiome_data.txt", "/home/docker/example_inputs/out.csv", "-s", "Sample", "-l", "Category", "-L", "3", "-n", "4", "--seed", "42", "--train_split", "0.8", "--compare_all_levels", "--nperm", "10") }
# these will be used by the argparser
opt <- load_args('taxaHFE-ML.R v1.0', 2)

## Run main ====================================================================

## set random seed
set_seed_func(opt$seed)

## check for outdir and make if not there
if (dir.exists(paste0(dirname(opt$OUTPUT))) == FALSE) {
  dir.create(path = paste0(dirname(opt$OUTPUT)))
}

## We need the flattened tree output, so always make it true
opt$write_flattened_tree = TRUE

## parameters specified
cat("\n","Parameters specified: \n\n")
cat("Global Options: \n")
cat(paste0("--subject_identifier: ", opt$subject_identifier), "\n")
cat(paste0("--label: ", opt$label), "\n")
cat(paste0("--feature_type: ", opt$feature_type), "\n")
cat(paste0("--cor_level: ", opt$cor_level), "\n")
cat(paste0("--ncores: ", opt$ncores), "\n")
cat(paste0("--seed: ", opt$seed), "\n")
cat("\nTaxaHFE Options: \n")
cat(paste0("--abundance: ", opt$abundance), "\n")
cat(paste0("--prevalence: ", opt$prevalence), "\n")
cat(paste0("--lowest_level: ", opt$lowest_level), "\n")
cat(paste0("--max_level: ", opt$max_level), "\n")
cat(paste0("--disable_super_filter: ", opt$disable_super_filter), "\n")
cat(paste0("--write_old_files: ", opt$write_old_files), "\n")
cat(paste0("--write_flattened_tree (always TRUE): ", opt$write_flattened_tree), "\n")
cat(paste0("--write_both_outputs: ", opt$write_both_outputs), "\n")
cat(paste0("--nperm: ", opt$nperm), "\n")
cat("\nDietML Options: \n")
cat(paste0("--train_split: ", opt$train_split), "\n")
cat(paste0("--model: ", opt$model), "\n")
cat(paste0("--folds: ", opt$folds), "\n")
cat(paste0("--metric: ", opt$metric), "\n")
cat(paste0("--tune_length: ", opt$tune_length), "\n")
cat(paste0("--tune_time: ", opt$tune_time), "\n")
cat(paste0("--tune_stop: ", opt$tune_stop), "\n")
cat(paste0("--shap: ", opt$shap), "\n")
cat(paste0("--compare-all-levels: ", opt$compare_all_levels), "\n")
cat(paste0("\n", "OUTPUT: ", opt$OUTPUT))

## check for inputs and read in read in ========================================
cat("\n\n", "###########################\n", "Reading in data...\n", "###########################")

## metadata file
metadata <- read_in_metadata(input = opt$METADATA, 
                             subject_identifier = opt$subject_identifier, 
                             label = opt$label)

## hierarchical data file ======================================================
hData <- read_in_hierarchical_data(input = opt$DATA, 
                                   metadata = metadata, 
                                   cores = opt$ncores)

## set initial test-train split ================================================
tr_te_split <- rsample::initial_split(metadata, prop = as.numeric(opt$train_split), strata = feature_of_interest)
train_metadata <- rsample::training(tr_te_split)
test_metadata  <- rsample::testing(tr_te_split)

##############################
## TAXAHFE FEATURE ENGINEERING
##############################

count = 1

for (split_metadata in list(train_metadata, test_metadata)) {
  
  ## Some messaging to let the user know we are performing TaxaHFE
  ## on Training AND Testing data
  if (count == 1) {
    cat("\n", "Generating Training Data:")
  } else if (count == 2) {
    cat("\n", "Generating Testing Data:")
  } else {
    stop("We have looped through your data too many times. This is weird. File a github issue.")
  }
  
  hData_split <- hData %>% dplyr::select(., clade_name, dplyr::any_of(split_metadata$subject_id))

  ## Build tree ================================================================
  cat("\n\n", "###########################\n", "Building Tree...\n", "###########################\n\n")
  cat("This may take a few minutes depending on how many features you have.\n")
  hTree <- build_tree(hData_split, 
                      filter_prevalence = opt$prevalence,
                      filter_mean_abundance = opt$abundance)

  ## Main competition ==========================================================
  cat("\n", "###########################\n", "Competing Tree...\n", "###########################\n\n")
  
  competed_tree <- compete_tree(
    hTree,
    # if in train loop (count = 1), compete to lowest level specified, else count = 2, barely compete, just one level
    lowest_level = switch(count, {opt$lowest_level}, { max(stringr::str_count(hData$clade_name, "\\|")) }),
    max_level = opt$max_level, # allows for all levels to be competed. Change to 1 for pairwise comparisons
    col_names = colnames(hData_split)[2:NCOL(hData_split)],
    # if in train loop (count = 1), corr competitions as specified, else count = 2, make almost everything a corr competition
    corr_threshold = switch(count, {opt$cor_level}, {as.numeric(0.1)}), 
    metadata = split_metadata,
    ncores = opt$ncores,
    feature_type = opt$feature_type,
    # if in train loop (count = 1), nperm as specified, else count = 2, barely permute the RF competitions
    nperm = switch(count, {opt$nperm}, {as.numeric(3)}), 
    disable_super_filter = opt$disable_super_filter
  )
  
  ## Extract information from tree  ============================================
  # Flatten the tree and tree decisions
  cat("\n", "############################################\n", "Flattening tree and writing final output...\n", "############################################\n\n")
  
  generate_outputs(
    competed_tree,
    split_metadata,
    colnames(hData_split)[2:NCOL(hData_split)],
    # if in train loop (count = 1), write opt$OUTPUT_train.csv, else count = 2, write opt$OUTPUT_test.csv
    switch(count, {paste0(tools::file_path_sans_ext(opt$OUTPUT), "_train.csv")}, {paste0(tools::file_path_sans_ext(opt$OUTPUT), "_test.csv")}), 
    opt$disable_super_filter,
    opt$write_both_outputs,
    opt$write_old_files,
    opt$write_flattened_tree,
    opt$ncores
  )
  
  ## iteratable to loop over train and test data
  count = count + 1
  
}

## message to user finish
cat(" Outputs written! TaxaHFE completed. \n")


##########################
## DIETML MACHINE LEARNING
##########################

if (opt$model == "none") {
  cat("Model set to none. DietML not run")
  save.image(file = paste0(dirname(opt$OUTPUT), "/taxaHFE_r_workspace.rds"))
} else { 
  
  ## format input for ML 
  if (opt$disable_super_filter) {
    train_data <- readr::read_csv(file = paste0(tools::file_path_sans_ext(opt$OUTPUT), "_train_no_sf.csv"))
  } else {
    train_data <- readr::read_csv(file = paste0(tools::file_path_sans_ext(opt$OUTPUT), "_train.csv"))
  }
  
  flattened_df_test <- readr::read_delim(file = paste0(tools::file_path_sans_ext(opt$OUTPUT), "_test_raw_data.tsv.gz"))
  
  ## make sure test and train have the same features
  flattened_df_test$name <- janitor::make_clean_names(flattened_df_test$name)
  test_data <- flattened_df_test %>%
    # only select rows with same feature names as train_data features
    dplyr::filter(., name %in% colnames(train_data)) %>%
    # only select columns in test_metadata
    dplyr::select(., name, dplyr::any_of(test_metadata$subject_id)) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(., var = "name") %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(., var = "subject_id")
  
  ## merge train and test data with metadata
  train_data <- merge((metadata %>% dplyr::select(., -feature_of_interest)), train_data, by = "subject_id")
  test_data <- merge(metadata, test_data, by = "subject_id")
  
  ## make colnames appropriate for ML (ranger is picky)
  colnames(train_data) <- make.names(colnames(train_data))
  colnames(test_data) <- make.names(colnames(test_data))
  
  ## organize the columns in test to be the same as train
  ## the manual test-train split (make_splits()) is picky
  ## they also have to have the same features
  ## The actual error that gets thrown (even if columns are just out of order):
  ## Error in `make_splits()` at scripts/models/dietML_ranger_tidy.R:41:1:
  ##   ! The analysis and assessment sets must have the same columns
  overlap_features <- dplyr::intersect(colnames(test_data), colnames(train_data))
  test_data <- test_data %>% dplyr::select(., dplyr::any_of(overlap_features))
  train_data_for_dietML <- train_data %>% dplyr::select(., dplyr::any_of(overlap_features))
  ## reorder test columns
  test_data_for_dietML <- test_data[names(train_data_for_dietML)]
  
  ## write the test and train data to file
  readr::write_csv(x = train_data_for_dietML, file = paste0(dirname(opt$OUTPUT), "/train_data.csv"))
  readr::write_csv(x = test_data_for_dietML, file = paste0(dirname(opt$OUTPUT), "/test_data.csv"))
  
  ## keep track of what type of analysis is being passed to ML
  opt$program <- "taxaHFE-ML"
  
  ## source dietML
  source("/scripts/dietML.R")

}

##########################
## COMPARE ALL LEVELS
##########################

if (opt$compare_all_levels) {
  opt$shap <- FALSE
  opt$model <- "rf"
  source("/scripts/utilities/compare_all_levels.R")
}
