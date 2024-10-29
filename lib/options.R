#!/usr/bin/env Rscript

## SCRIPT: options.R ===================================================
## DATE:   June, 24 2024
##
## PURPOSE: Holds the commandline args for taxaHFE
## the options.R script

## docker info =================================================================
suppressPackageStartupMessages(library(argparse, quietly = T, verbose = F, warn.conflicts = F))

# Function to add base arguments common to all programs
add_taxa_hfe_base_args <- function(parser) {
  taxa_hfe_base_args = parser$add_argument_group('TaxaHFE arguments', 'Options to pass to TaxaHFE')
  taxa_hfe_base_args$add_argument('-s', '--subject_identifier', type="character", metavar="<string>", default='subject_id', help='Metadata column name containing subject IDs')
  taxa_hfe_base_args$add_argument('-l', '--label', type="character", metavar="<string>", default='feature_of_interest', help='Metadata column name of interest for ML')
  taxa_hfe_base_args$add_argument('-t', '--feature_type', type="character", metavar="<string>", default='factor', help='Is the ML label a factor or numeric')
  taxa_hfe_base_args$add_argument('-a', '--abundance', type="numeric", metavar="<numeric>", default='0', help='Minimum mean abundance of feature')
  taxa_hfe_base_args$add_argument('-p', '--prevalence', type="numeric", metavar="<numeric>", default='0.01', help='Minimum prevalence of feature')
  taxa_hfe_base_args$add_argument('-L', '--lowest_level', type="integer", metavar="<numeric>", default='2', help='Most general level allowed to compete')
  taxa_hfe_base_args$add_argument('-m', '--max_level', type="integer", metavar="<numeric>", default='1000', help='How many hierarchical levels should be allowed to compete')
  taxa_hfe_base_args$add_argument('-c', '--cor_level', type="numeric", metavar="<numeric>", default='0.95', help='Initial pearson correlation filter')
  taxa_hfe_base_args$add_argument('-d', '--disable_super_filter', action="store_true", help='Disable running of the super filter (final forest competition)')
  taxa_hfe_base_args$add_argument('-w', '--write_old_files', action="store_true", help='Write individual level files and old HFE files')
  taxa_hfe_base_args$add_argument('-W', '--write_flattened_tree', action="store_true", help='Write a compressed backup of the entire competed tree')
  taxa_hfe_base_args$add_argument('-D', '--write_both_outputs', action="store_true", help='Write an output for pre and post super filter results, overridden by --disable_super_filter')
  taxa_hfe_base_args$add_argument('--nperm', type="integer", metavar="<numeric>", default='40', help='Number of taxaHFE RF permutations')
  parser$add_argument('-n', '--ncores', type="integer", metavar="<numeric>", default='2', help='Number of CPU cores to use')
  parser$add_argument('--seed', type="numeric", metavar="<numeric>", help='Set a random numeric seed. If None, defaults to system time')
  return(parser)
}

# Function to initialize parser for a program
initialize_parser <- function(version, program_name, description) {
  parser <- argparse::ArgumentParser(
    description=description,
    usage=paste(program_name, "[options] METADATA DATA OUTPUT"),
    formatter_class="type('CustomFormatter', (argparse.ArgumentDefaultsHelpFormatter, argparse.MetavarTypeHelpFormatter, argparse.RawTextHelpFormatter), {})"
  )
  
  # Common arguments for all programs
  parser$add_argument('METADATA', metavar='METADATA', type="character", help="path to metadata input (txt | tsv | csv)")
  parser$add_argument('DATA', metavar='DATA', type="character", help="path to input file from hierarchical data (i.e. hData data) (txt | tsv | csv)")
  parser$add_argument('OUTPUT', metavar='OUTPUT', type="character", help="output file name (csv)")
  
  parser$add_argument('-v', '--version', action='version', version=version)
  parser$add_argument('--data_dir', type="character", metavar="<string>", default='/data', help='Directory for data files input/output, ignored if using absolute paths')
  
  return(parser)
}

# Function to add taxaHFE-ML specific arguments
add_taxa_hfe_ml_args <- function(parser) {
  taxa_hfe_ml_base_args = parser$add_argument_group('TaxaHFE-ML specific arguments', 'Options to pass to TaxaHFE-ML for machine learning and SHAP analysis of TaxaHFE features')
  taxa_hfe_ml_base_args$add_argument('--train_split', type="numeric", metavar="<numeric>", default='0.8', help='Percentage of samples to use for training')
  taxa_hfe_ml_base_args$add_argument('--model', type="character", metavar="<string>", default='rf', choices=c("rf", "enet", "none"), help='ML model to use')
  taxa_hfe_ml_base_args$add_argument('--folds', type="numeric", metavar="<numeric>", default='10', help='Number of CV folds for tuning')
  taxa_hfe_ml_base_args$add_argument('--metric', type="character", metavar="<string>", default='bal_accuracy', choices=c("roc_auc", "bal_accuracy", "accuracy", "mae", "rmse", "rsq", "kap", "f_meas", "ccc"), help='Metric to optimize')
  taxa_hfe_ml_base_args$add_argument('--tune_length', type="numeric", metavar="<numeric>", default='80', help='Number of hyperparameter combinations to sample')
  taxa_hfe_ml_base_args$add_argument('--tune_time', type="numeric", metavar="<numeric>", default='2', help='Time for hyperparameter search (in hours)')
  taxa_hfe_ml_base_args$add_argument('--tune_stop', type="numeric", metavar="<numeric>", default='10', help='Number of HP iterations without improvement before stopping')
  taxa_hfe_ml_base_args$add_argument('--permute', type="numeric", metavar="<numeric>", default='1', help='Number of times to permute the ML assessment process, resulting in n different test/train split inputs')
  taxa_hfe_ml_base_args$add_argument('--shap', action="store_true", help='Calculate SHAP values')
  return(parser)
}

# Function to add taxaHFE-RM specific arguments
add_taxa_hfe_rm_args <- function(parser) {
  taxa_hfe_rm_base_args = parser$add_argument_group('TaxaHFE-RM specific arguments', 'Options to pass to DietML for machine learning and SHAP analysis of TaxaHFE features')
  #taxa_hfe_rm_base_args$add_argument('--train_split', type="numeric", metavar="<numeric>", default='0.8', help='Percentage of samples to use for training')
  #taxa_hfe_rm_base_args$add_argument('--model', type="character", metavar="<string>", default='rf', choices=c("rf", "enet", "none"), help='ML model to use')
  #taxa_hfe_rm_base_args$add_argument('--shap', action="store_true", help='Calculate SHAP values')

  return(parser)
}

## load args function
load_args <- function(program) {
  version <- Sys.getenv('TAXA_HFE_VERSION')
  if (version == "") {
    version <- "0"
  }

  if (program == "taxa_hfe") {
    # Initialize the parser for taxaHFE
    parser <- initialize_parser(version, program, "Hierarchical feature engineering (HFE) for feature reduction")
    parser <- add_taxa_hfe_base_args(parser)
  } else if (program == "taxa_hfe_ml") {
    # Initialize the parser for taxaHFE-ML
    parser <- initialize_parser(version, program, "Hierarchical feature engineering (HFE) with ML")
    parser <- add_taxa_hfe_base_args(parser)
    parser <- add_taxa_hfe_ml_args(parser)
  } else if (program == "taxa_hfe_rm") {
    # Initialize the parser for taxaHFE-RM (you can define specific args for this if needed)
    parser <- initialize_parser(version, program, "Hierarchical feature engineering (HFE) with RM")
    parser <- add_taxa_hfe_base_args(parser)
    # Add any RM-specific arguments here if necessary, e.g.:
    # parser <- add_taxa_hfe_rm_args(parser)
  } else {
    stop("Invalid program name. Choose from 'taxa_hfe', 'taxa_hfe_ml', or 'taxa_hfe_rm'.")
  }
  
  # Parse the command-line arguments and return them
  opt <- parser$parse_args(commandArgs(TRUE))

  # also normalize all input paths to the data_dir
  # will ignore the data_dir if the path is abosolute
  opt$METADATA <- file.path(opt$data_dir, opt$METADATA)
  opt$DATA <- file.path(opt$data_dir, opt$DATA)
  opt$OUTPUT <- file.path(opt$data_dir, opt$OUTPUT)

  return(opt)
}
