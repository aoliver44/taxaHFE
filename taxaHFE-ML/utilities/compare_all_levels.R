#!/usr/bin/env Rscript

## SCRIPT: compare_all_levels.R ===================================================
## AUTHOR: Andrew Oliver
## DATE:   Jun, 18 2024
##
## PURPOSE: Compare RF metrics taxaHFE-ML, taxaHFE(+/- SF) and all summarized levels
## the compare_all_levels.R script

## docker info =================================================================
#docker run --rm -v `pwd`:/home/docker -w /home/docker aoliver44/leakage_free_taxaHFE:latest

cat("\n", "############################################\n", "Comparing all levels...\n", "############################################\n\n")

## Run taxaHFE(+/- SF) first

## Build tree ==================================================================
hTree <- build_tree(hData, 
                    filter_prevalence = opt$prevalence,
                    filter_mean_abundance = opt$abundance)

## Main competition ============================================================
competed_tree <- compete_tree(
  hTree,
  lowest_level = opt$lowest_level,
  max_level = opt$max_level, # allows for all levels to be competed. Change to 1 for pairwise comparisons
  col_names = colnames(hData)[2:NCOL(hData)],
  corr_threshold = opt$cor_level,
  metadata = metadata,
  ncores = opt$ncores,
  feature_type = opt$feature_type,
  nperm = opt$nperm,
  disable_super_filter = FALSE
)

## Extract information from tree  ==============================================
# Flatten the tree and tree decisions
generate_outputs(
  competed_tree,
  metadata,
  colnames(hData)[2:NCOL(hData)],
  opt$OUTPUT, opt$disable_super_filter,
  write_both_outputs = TRUE,
  write_old_files = TRUE,
  write_flattened_df_backup = FALSE,
  opt$ncores
)

max_levels <- max(stringr::str_count(hData$clade_name, "\\|")) + 1
files <- c(paste0(tools::file_path_sans_ext(opt$OUTPUT), ".csv"), 
           paste0(tools::file_path_sans_ext(opt$OUTPUT), "_no_sf.csv"), 
           paste0(paste0(tools::file_path_sans_ext(opt$OUTPUT),"_level_", seq(1:max_levels), ".csv")))
program <- c("taxaHFE_SF", "taxaHFE", paste0("level_",seq(1:max_levels)))
count = 1

for (file in files) {
  
  cat(paste0("Competeting level: "), opt$program)
  ## read in raw data
  raw_compare_df <- readr::read_csv(files[count])
  opt$program <- program[count]
  
  ## split data as taxaHFE-ML was split
  train_data_for_dietML <- raw_compare_df %>% dplyr::filter(., subject_id %in% train_metadata$subject_id)
  test_data_for_dietML <- raw_compare_df %>% dplyr::filter(., subject_id %in% test_metadata$subject_id)
  
  source("/scripts/dietML.R")
  
  count = count + 1
}




