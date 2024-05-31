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
# ex. commandArgs <- function(x) { c("example_inputs/metadata.txt", "example_inputs/microbiome_data.txt", "example_inputs/out.txt", "-s", "Sample", "-l", "Category", "-L", "3", "-n", "4", "--seed", "42", "--shap", "--train_split", "0.8") }
# these will be used by the argparser
opt <- load_args('taxaHFE-ML.R v1.0', 2)
## Run main ====================================================================

## set random seed
set_seed_func(opt$seed)

## check for outdir and make if not there
if (dir.exists(paste0(dirname(opt$OUTPUT))) == FALSE) {
  dir.create(path = paste0(dirname(opt$OUTPUT)))
}

## We need the flattened tree outout, so always make it true
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
cat(paste0("--max_depth: ", opt$max_depth), "\n")
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
    lowest_level = switch(count, {opt$lowest_level}, {max(flattened_df_train$depth) - 1}),
    max_depth = opt$max_depth, # allows for all levels to be competed. Change to 1 for pairwise comparisons
    col_names = colnames(hData_split)[2:NCOL(hData_split)],
    corr_threshold = switch(count, {opt$cor_level}, {as.numeric(0.1)}), # massively reduce the RF compititons by decreasing corr_thresh - run only on test data - speed up
    metadata = split_metadata,
    ncores = opt$ncores,
    feature_type = opt$feature_type,
    nperm = switch(count, {opt$nperm}, {as.numeric(3)}), # massively reduce the RF nperm by decreasing nperm - run only on test data - speed up
    disable_super_filter = opt$disable_super_filter
  )
  
  ## Extract information from tree  ============================================
  # Flatten the tree and tree decisions
  cat("\n", "############################################\n", "Flattening tree and writing final output...\n", "############################################\n\n")
  
  generate_outputs(
    competed_tree,
    split_metadata,
    colnames(hData_split)[2:NCOL(hData_split)],
    opt$OUTPUT, opt$disable_super_filter,
    opt$write_both_outputs,
    opt$write_old_files,
    opt$write_flattened_tree,
    opt$ncores
  )
  
  ## rename output objects for ease
  ## if user disabled superfilter, our outputs will reflect that
  if (count == 1) {
    assign(x = "flattened_df_train", value = flattened_df, envir = .GlobalEnv)
    if (opt$disable_super_filter == TRUE) {
      assign(x = "train_data", output_1_no_sf, envir = .GlobalEnv)
    } else {
      assign(x = "train_data", output_1, envir = .GlobalEnv)
    }
  } else {
    assign(x = "flattened_df_test", value = flattened_df, envir = .GlobalEnv)
    if (opt$disable_super_filter == TRUE) {
      assign(x = "test_data", output_2_no_sf, envir = .GlobalEnv)
    } else {
      assign(x = "test_data", output_2, envir = .GlobalEnv)
    }
  }
  
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
  source("/scripts/dietML.R")
}