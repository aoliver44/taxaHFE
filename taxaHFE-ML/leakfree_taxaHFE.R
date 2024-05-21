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

## source scripts ==============================================================
source("/scripts/tree.R")

## add commandline options =====================================================

'Hierarchical feature engineering (HFE) for the reduction of features with respects to a factor or regressor
Usage:
    taxaHFE-ML [options] <METADATA> <DATA> <OUTPUT>

Global Options:
    -h --help                         Show help text.
    -v --version                      Show version.
    -s --subject_identifier <string>  Metadata column name containing subject IDs [default: subject_id]
    -l --label <string>               Metadata column name of interest for ML [default: cluster]
    -t --feature_type <string>        Is the ML label a factor or numeric [default: factor]
    -c --cor_level <float>            Initial pearson correlation filter [default: 0.95]
    -n --ncores <int>                 Number of cpu cores to use [default: 2]
    --seed <numeric>                  Set a random numeric seed, default is to use system time
TaxaHFE Options:
    -a --abundance <float>            Minimum mean abundance of feature [default: 0.0001]
    -p --prevalence <float>           Minimum prevalence of feature [default: 0.01]
    -L --lowest_level <int>           Most general level allowed to compete [default: 2]
    -m --max_depth <int>              How many hierarchical levels should be allowed to compete [default: 1000]
    -d --disable_super_filter         Disable running of the super filter (final forest competition)
    -w --write_old_files              Write individual level files and old HFE files
    -W --write_flattened_tree         Write a compressed backup of the entire competed tree
    -D --write_both_outputs           Write an output for pre and post super filter results, overridden by --disable_super_filter
    --nperm <int>                     Number of RF permutations [default: 40]
DietML Options:
    --train_split what percentage of samples should be used in training 
            [default: 0.70]
    --model what model would you like run 
            (options: rf,enet) [default: rf]
    --folds number of CV folds to tune with [default: 10]
    --metric what metric would you like to optimize in training 
            (options: roc_auc, bal_accuracy, accuracy, mae, rmse, rsq, kap, 
             f_meas, ccc) [default: bal_accuracy]
    --tune_length number of hyperparameter combinations to sample [default: 80]
    --tune_time length of time tune_bayes runs [default: 10]
    --tune_stop number of HP interations to let pass without a metric 
            improvement [default: 10]
    --shap attempt to calcualte shap values? [default: TRUE]

Arguments:
    METADATA path to metadata input (txt | tsv | csv)
    DATA path to input file from hierarchical data (i.e. hData data) (txt | tsv | csv)
    OUTPUT output file name (csv)

' -> doc

# these options will be converted to numeric by load_docopt
numeric_options <- c("train_split", "abundance", "prevalence", "lowest_level", "max_depth", "cor_level", "ncores", "nperm", "folds", "tune_length", "tune_time", "tune_stop", 
                     "--train_split", "--abundance", "--prevalence", "--lowest_level", "--max_depth", "--cor_level", "--ncores", "--nperm", "--folds", "--tune_length", "--tune_time", "--tune_stop")
# to use this code line-by-line in the Rstudio context, commandArgs can be overloaded to specify the desired flags
# ex. commandArgs <- function(x) { "-s subject_id -l new_butyrate -L 3 -n 4 -wWD --seed 42 /home/docker/metadata.txt example_inputs/microbiome_data.txt example_inputs/out.csv" }
# these will be used by the options loader
#commandArgs <- function(x) { "-s subject_id -l crp_boxcox -t numeric -n 4 -a 0 -L 3 -d --train_split 0.7 --seed 56 --metric mae /home/docker/example_data/stephanie_crp_metadata.csv /home/docker/example_data/input_file_relabund.txt /home/docker/ml_output/output.csv" }
opt <- load_docopt(doc, version = 'taxaHFE-ML.R v1.0\n\n', to_convert = numeric_options)

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