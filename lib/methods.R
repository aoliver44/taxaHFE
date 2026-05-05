#!/usr/bin/env Rscript

## SCRIPT: methods.R ===================================================
## DATE:   June, 24 2024
##
## PURPOSE: Holds the methods for taxaHFE and associated programs
## the methods.R script

## docker info =================================================================

method_taxa_hfe <- function(
  h_data, metadata, prevalence, abundance,
  lowest_level, max_level, cor_level, ncores,
  feature_type, nperm, disable_super_filter,
  write_both_outputs, write_flattened_tree, write_old_files,
  col_names, output, random_effects
) {
  ## Build tree ================================================================
  cat("\n\n", "###########################\n", "Building Tree...\n", "###########################\n\n")
  cat("This may take a few minutes depending on how many features you have.\n")
  h_tree <- build_tree(
    h_data,
    filter_prevalence = prevalence,
    filter_mean_abundance = abundance
  )

  ## Main competition ==========================================================
  competed_tree <- compete_tree(
    tree = h_tree,
    lowest_level = lowest_level,
    max_level = max_level, # allows for all levels to be competed. Change to 1 for pairwise comparisons
    col_names = col_names,
    corr_threshold = cor_level,
    metadata = metadata,
    ncores = ncores,
    feature_type = feature_type,
    nperm = nperm,
    disable_super_filter = disable_super_filter, 
    random_effects = random_effects
  )

  ## write outputs =============================================================
  generate_outputs(
    tree = competed_tree,
    metadata = metadata,
    col_names = col_names,
    output_location = output,
    disable_super_filter = disable_super_filter,
    write_both_outputs = write_both_outputs,
    write_old_files = write_old_files,
    write_flattened_df_backup = write_flattened_tree,
    ncores = ncores
  )
}

method_taxa_hfe_ml <- function(
  h_data, metadata, prevalence, abundance,
  lowest_level, max_level, cor_level, ncores,
  feature_type, nperm, disable_super_filter, train_metadata,
  test_metadata, seed, random_effects, write_flattened_tree, 
  output
) {
  ## TRAIN DATA ================================================================
  ## training data - gets HFE treatment
  ## subset h_data for training samples
  training_h_data_split <- h_data %>% dplyr::select(., clade_name, dplyr::any_of(train_metadata$subject_id))
  ## build training tree
  h_tree_train <- build_tree(training_h_data_split, filter_prevalence = prevalence, filter_mean_abundance = abundance)
  ## compete training tree
  competed_tree <- compete_tree(
    h_tree_train,
    lowest_level = lowest_level,
    max_level = max_level, 
    col_names = colnames(training_h_data_split)[2:NCOL(training_h_data_split)],
    corr_threshold = cor_level,
    metadata = train_metadata,
    ncores = ncores,
    feature_type = feature_type,
    nperm = nperm,
    disable_super_filter = disable_super_filter,
    random_effects = random_effects
  )
  ## Extract information from tree  ============================================
  train_data <- prepare_flattened_df(
    train = TRUE,
    levels = FALSE,
    node = competed_tree,
    metadata = train_metadata,
    disable_super_filter = disable_super_filter,
    col_names = colnames(training_h_data_split)[2:NCOL(training_h_data_split)],
    write_flattened_tree = write_flattened_tree,
    ncores = ncores,
    output = output
  )
  ## TEST DATA =================================================================
  ## subset h_data for testing samples 
  testing_h_data_split <- h_data %>% dplyr::select(., clade_name, dplyr::any_of(test_metadata$subject_id))
  ## build testing tree - no filters get run! ALL testing features make it through
  h_tree_test <- build_tree(testing_h_data_split, filter_prevalence = 0, filter_mean_abundance = 0)
  ## Extract information from tree 
  test_data <- prepare_flattened_df(
    train = FALSE,
    levels = FALSE,
    node = h_tree_test,
    metadata = test_metadata,
    disable_super_filter = TRUE, ## same as abundance and prev filters, ALL testing features make it through
    col_names = colnames(testing_h_data_split)[2:NCOL(testing_h_data_split)],
    write_flattened_tree = FALSE,
    output = output
  )

  ## make colnames appropriate for ML (ranger is picky)
  colnames(train_data) <- make.names(colnames(train_data))
  colnames(test_data) <- make.names(colnames(test_data))

  ## organize the columns in test to be the same as train
  ## the manual test-train split (make_splits()) is picky
  ## they also have to have the same features
  ## The actual error that gets thrown (even if columns are just out of order):
  ## Error in `make_splits()` at lib/models/diet_ml_ranger_tidy.R:41:1:
  ##   ! The analysis and assessment sets must have the same columns
  test_data <- test_data %>% dplyr::select(., dplyr::any_of(colnames(train_data)))
  ## reorder test columns
  test_data <- test_data[names(train_data)]
  
  ## Error out if training data columns and testing columns do not match
  ## (if working well, should be TRUE. Below it is asking if it is the inverse
  ## of it is TRUE...basically a problem)
  if (!all(colnames(train_data)==colnames(test_data))) {
    logger::log_fatal("The training colnames and testing colnames do not match. This is unusual -
         take note of the command and your data and raise an issue on our GitHub
         (https://github.com/aoliver44/taxaHFE/issues)")
    stop()
  }

  ## store data in list of data for dietML
  diet_ml_inputs <- list(train_data, test_data)

  return(diet_ml_inputs)
}
  

method_levels <- function(
  h_data, metadata, ncores, 
  prevalence, abundance,
  disable_super_filter,
  seed
) {

  ## Build tree ================================================================
  cat("\n\n", "###########################\n", "Building levels tree...\n", "###########################\n\n")
  h_tree_levels <- build_tree(
    h_data,
    filter_prevalence = prevalence,
    filter_mean_abundance = abundance
  )

  ## flatten the data
  flattened_df <- prepare_flattened_df(
    train = FALSE,
    levels = TRUE,
    node = h_tree_levels,
    metadata = metadata,
    disable_super_filter = disable_super_filter,
    col_names = colnames(h_data)[2:NCOL(h_data)],
    write_flattened_tree = FALSE,
    ncores = ncores,
    output = output
  )

  ## attach the summary files to dietML_input list
  diet_ml_level_inputs <- generate_summary_files(
    input = flattened_df,
    metadata = metadata,
    disable_super_filter = ifelse(disable_super_filter, "no_sf", "sf"),
    seed = seed
  )

  return(diet_ml_level_inputs)
}
