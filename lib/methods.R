#!/usr/bin/env Rscript

## SCRIPT: methods.R ===================================================
## DATE:   June, 24 2024
##
## PURPOSE: Holds the methods for taxaHFE and associated programs
## the methods.R script

## docker info =================================================================

method_taxa_hfe <- function(hdata, metadata, prevalence, abundance,
                           lowest_level, max_level, cor_level, ncores,
                           feature_type, nperm, disable_super_filter,
                           write_both_outputs, write_flattened_tree, col_names,
                           target_list, output, seed, random_effects) {
  ## Build tree ================================================================
  cat("\n\n", "###########################\n", "Building Tree...\n", "###########################\n\n")
  cat("This may take a few minutes depending on how many features you have.\n")
  hTree <- build_tree(hData,
    filter_prevalence = prevalence,
    filter_mean_abundance = abundance
  )

  ## Main competition ==========================================================
  competed_tree <- compete_tree(
    hTree,
    lowest_level = lowest_level,
    max_level = max_level, # allows for all levels to be competed. Change to 1 for pairwise comparisons
    col_names = colnames(hData)[2:NCOL(hData)],
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
    competed_tree,
    metadata,
    colnames(hData)[2:NCOL(hData)],
    opts$OUTPUT, opts$disable_super_filter,
    opts$write_both_outputs,
    opts$write_old_files,
    opts$write_flattened_tree,
    opts$ncores
  )
  
  ## Extract information from tree  ============================================
  ## store data in list of data for dietML
  
  flattened_df <- prepare_flattened_df(node = competed_tree,
                                       metadata = metadata,
                                       disable_super_filter = disable_super_filter,
                                       col_names = colnames(hData)[2:NCOL(hData)]
  )
  
  ## store taxaHFE outputs in list
  diet_ml_inputs <<- store_diet_ml_inputs(target_list = diet_ml_inputs,
    object = flattened_df,
    super_filter = ifelse(disable_super_filter, "no_sf", "sf"),
    method = "taxa_hfe",
    train_test_attr = NA,
    level_n = NA,
    seed = seed
    )
  
}

method_taxa_hfe_ml <- function(hdata, metadata, prevalence, abundance,
                              lowest_level, max_level, cor_level, ncores,
                              feature_type, nperm, disable_super_filter,
                              write_both_outputs, write_flattened_tree,
                              train_split, model, folds, metric, tune_length,
                              tune_time, tune_stop, shap, target_list, output, seed, random_effects) {

  count <- 1

  for (split_metadata in list(train_metadata, test_metadata)) {
    ## Some messaging to let the user know we are performing TaxaHFE
    hData_split <- hData %>% dplyr::select(., clade_name, dplyr::any_of(split_metadata$subject_id))

    ## Build tree ================================================================
    hTree <- build_tree(hData_split,
                        filter_prevalence = prevalence,
                        filter_mean_abundance = abundance
    )

    ## Main competition ==========================================================
    competed_tree <- compete_tree(
      hTree,
      # if in train loop (count = 1), compete to lowest level specified, else count = 2, barely compete, just one level
      lowest_level = switch(count, {lowest_level}, {max(stringr::str_count(hData$clade_name, "\\|"))}),
      max_level = max_level, # allows for all levels to be competed. Change to 1 for pairwise comparisons
      col_names = colnames(hData_split)[2:NCOL(hData_split)],
      # if in train loop (count = 1), corr competitions as specified, else count = 2, make almost everything a corr competition
      corr_threshold = switch(count, {cor_level}, {as.numeric(0.1)}),
      metadata = split_metadata,
      ncores = ncores,
      feature_type = feature_type,
      # if in train loop (count = 1), nperm as specified, else count = 2, barely permute the RF competitions
      nperm = switch(count, {nperm}, {as.numeric(3)}),
      disable_super_filter = disable_super_filter,
      random_effects = random_effects
    )

    ## Extract information from tree  ============================================
    if (count == 1) {
      train_data <- prepare_flattened_df(node = competed_tree,
        metadata = metadata,
        disable_super_filter = disable_super_filter,
        col_names = colnames(hData_split)[2:NCOL(hData_split)])
    } else {
      flattened_df_test <- flatten_tree_with_metadata(competed_tree)
      col_names = colnames(hData_split)[2:NCOL(hData_split)]
      colnames(flattened_df_test)[11:NCOL(flattened_df_test)] <- col_names
      
      ## clean pathString names and use these, they will always be unique
      ## downside, longer names
      flattened_df_test$pathString <- flattened_df_test$pathString %>% janitor::make_clean_names()
      flattened_df_test$pathString <- gsub(pattern = "taxa_tree_", replacement = "", x = flattened_df_test$pathString)
      
      test_data <- flattened_df_test %>%
        dplyr::select(., pathString, 11:dplyr::last_col()) %>%
        tibble::remove_rownames() %>%
        tibble::column_to_rownames(., var = "pathString") %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "subject_id")
      
      test_data <- merge(metadata, test_data, by = "subject_id")
    }
    
  ## iteratable to loop over train and test data
  count <- count + 1
  }

  ## make colnames appropriate for ML (ranger is picky)
  colnames(train_data) <- make.names(colnames(train_data))
  colnames(test_data) <- make.names(colnames(test_data))

  ## organize the columns in test to be the same as train
  ## the manual test-train split (make_splits()) is picky
  ## they also have to have the same features
  ## The actual error that gets thrown (even if columns are just out of order):
  ## Error in `make_splits()` at lib/models/diet_ml_ranger_tidy.R:41:1:
  ##   ! The analysis and assessment sets must have the same columns
  overlap_features <- dplyr::intersect(colnames(test_data), colnames(train_data))
  test_data <- test_data %>% dplyr::select(., dplyr::any_of(overlap_features))
  train_data_for_diet_ml <- train_data %>% dplyr::select(., dplyr::any_of(overlap_features))
  ## reorder test columns
  test_data_for_diet_ml <- test_data[names(train_data_for_diet_ml)]

  ## store data in list of data for dietML
  diet_ml_inputs <<- store_diet_ml_inputs(target_list = diet_ml_inputs,
    object = train_data_for_diet_ml,
    super_filter = ifelse(disable_super_filter, "no_sf", "sf"),
    method = "taxa_hfe_ml",
    train_test_attr = "train",
    level_n = NA,
    seed = seed)
  diet_ml_inputs <<- store_diet_ml_inputs(target_list = diet_ml_inputs,
    object = test_data_for_diet_ml,
    super_filter = ifelse(disable_super_filter, "no_sf", "sf"),
    method = "taxa_hfe_ml",
    train_test_attr = "test",
    level_n = NA,
    seed = seed)

}
  

method_levels <- function(hdata, metadata, prevalence, abundance,
                          lowest_level, max_level, cor_level, ncores,
                          feature_type, nperm, disable_super_filter,
                          write_both_outputs, write_flattened_tree, col_names,
                          target_list, output, seed, random_effects) {
  
  ## Build tree ================================================================
  cat("\n\n", "###########################\n", "Building Tree...\n", "###########################\n\n")
  cat("This may take a few minutes depending on how many features you have.\n")
  hTree <- build_tree(hData,
                      filter_prevalence = prevalence,
                      filter_mean_abundance = abundance
  )
  
  ## Main competition ==========================================================
  ## fixed some competition parameters to make it much faster (for the levels,
  ## really the tree building part is important, not the actual competition)
  competed_tree <- compete_tree(
    hTree,
    lowest_level = max(stringr::str_count(hData$clade_name, "\\|")),
    max_level = max_level, # allows for all levels to be competed. Change to 1 for pairwise comparisons
    col_names = colnames(hData)[2:NCOL(hData)],
    corr_threshold = cor_level,
    metadata = metadata,
    ncores = ncores,
    feature_type = feature_type,
    nperm = 3, # hard encoded because the levels dont need an RF competition
    disable_super_filter = disable_super_filter,
    random_effects = random_effects
  )
  
  ## flatten the data
  flattened_df <- flatten_tree_with_metadata(competed_tree)
  ## add back in the col names
  colnames(flattened_df)[11:NCOL(flattened_df)] <- col_names
  
  ## attach the summary files to dietML_input list
  generate_summary_files(input = flattened_df, 
                         metadata = metadata, 
                         target_list = diet_ml_inputs, 
                         disable_super_filter = ifelse(disable_super_filter, "no_sf", "sf"),
                         seed = seed
                        )
}

