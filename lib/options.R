source("lib/validators.R")

suppressPackageStartupMessages(library(argparse, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(parallel, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE))

# generates a seed for the default value of --seed
# seeds must be valid R 32-bit integers, so the range is between (-max int, max int)
default_seed <- function() {
  return(as.integer(runif(1, min = -1 * .Machine$integer.max, max = .Machine$integer.max)))
}

# these are the argument groups, each list corresponds to an argument grouping
# name and description are used to generate the argument group
# args represents the group of actual flags, the function arguments passed to parser$add_argument
# each of these is grouped to be used inside a parser$add_argument_group
# these are pulled outside into objects so they can be referenced during testing and easily added to
# to add a new parse group
# 1. Add a new group to this list, including a name, desc, and group of actual args
# 2. Add a new method for loading the args, ex. load_taxa_hfe_ml_args, all these do is hold the desired argument groups and call return(load_args(...))
# 3. Add a new entry to validators for any flag that needs validation
# 4. Add a new item to test_flag_values in test_options.R, including an non-default value for each flag to test, and any errors/warnings that should be tested as well
# 5. Add a new item to parser_flag_values_map in test_options.R that references the new load_taxa_hfe_*** method and the desired flag values
argument_groups <- list(
  taxa_hfe_base_args=list(
    name="TaxaHFE arguments",
    desc="Options to pass to TaxaHFE",
    args=list(
      subject_identifier=list("-s", "--subject_identifier", type="character", metavar="<string>", default="subject_id", help="Metadata column name containing subject IDs"),
      label=list("-l", "--label", type="character", metavar="<string>", default="feature_of_interest", help="Metadata column name of interest for ML"),
      feature_type=list("-t", "--feature_type", type="character", metavar="<string>", default="factor", help="Is the ML label a factor or numeric"),
      random_effects=list("-R", "--random_effects", action="store_true", help="Consider repeated measures. Note: columns 'individual' and 'time' must be present. [BETA]"),
      k_splits=list("-k", "--k_splits", type="numeric", metavar="<numeric>", default="3", help="We use kmeans to factorize a numeric response for repeated measures. How many categories should we create? [BETA]"),
      abundance=list("-a", "--abundance", type="numeric", metavar="<numeric>", default="0", help="Minimum mean abundance of feature"),
      prevalence=list("-p", "--prevalence", type="numeric", metavar="<numeric>", default="0.01", help="Minimum prevalence of feature"),
      lowest_level=list("-L", "--lowest_level", type="integer", metavar="<numeric>", default="3", help="Most general level allowed to compete"),
      max_level=list("-m", "--max_level", type="integer", metavar="<numeric>", default="15", help="How many hierarchical levels should be allowed to compete"),
      cor_level=list("-c", "--cor_level", type="numeric", metavar="<numeric>", default="0.95", help="Initial pearson correlation filter"),
      disable_super_filter=list("-d", "--disable_super_filter", action="store_true", help="Disable running of the super filter (final forest competition)"),
      write_old_files=list("-w", "--write_old_files", action="store_true", help="Write individual level files and old HFE files"),
      write_flattened_tree=list("-W", "--write_flattened_tree", action="store_true", help="Write a compressed backup of the entire competed tree"),
      write_both_outputs=list("-D", "--write_both_outputs", action="store_true", help="Write an output for pre and post super filter results, overridden by --disable_super_filter"),
      nperm=list("--nperm", type="integer", metavar="<numeric>", default="40", help="Number of taxaHFE RF permutations"),
      ncores=list("-n", "--ncores", type="integer", metavar="<numeric>", default="2", help="Number of parallel processes to run in certain portions of taxaHFE that support parallel processing. To limit overall resource usage of taxaHFE, limit the amount of resources available to the container (e.g. --cpus=4 for Docker)")
    )
  ),
  taxa_hfe_ml_args=list(
    name="TaxaHFE-ML specific arguments",
    desc="Options to pass to TaxaHFE-ML for machine learning and SHAP analysis of TaxaHFE features",
    args=list(
      train_split=list("--train_split", type="numeric", metavar="<numeric>", default="0.8", help="Percentage of samples to use for training"),
      model=list("--model", type="character", metavar="<string>", default="rf", choices=c("rf", "enet"), help="ML model to use"),
      folds=list("--folds", type="numeric", metavar="<numeric>", default="10", help="Number of CV folds for tuning"),
      metric=list("--metric", type="character", metavar="<string>", default="bal_accuracy", choices=c("roc_auc", "bal_accuracy", "accuracy", "mae", "rmse", "rsq", "kap", "f_meas", "ccc"), help="Metric to optimize"),
      tune_length=list("--tune_length", type="numeric", metavar="<numeric>", default="80", help="Number of hyperparameter combinations to sample"),
      tune_time=list("--tune_time", type="numeric", metavar="<numeric>", default="2", help="Time for hyperparameter search (in minutes)"),
      tune_stop=list("--tune_stop", type="numeric", metavar="<numeric>", default="10", help="Number of HP iterations without improvement before stopping"),
      permute=list("--permute", type="numeric", metavar="<numeric>", default="1", help="Number of times to permute the ML assessment process, resulting in n different test/train split inputs"),
      shap=list("--shap", action="store_true", help="Calculate SHAP values"),
      summarized_levels=list("--summarized_levels", action="store_true", help="Include summarized levels in ML competition")
    )
  )
)

# flag name mapped to a validator function
# validator functions always have the signature function(flag_name, flag_value, opts)
# this method should stop(error_message) when they are not satisfied
# or warning(message)
# ex.
# validate_some_flag <- function(flag_name, flag_value, opts) {
#   # check values as needed
#   # all_flags is fully parsed 
#   if (error) {
#     stop(error_message)
#   } else if (warning) {
#     warning(warning_message)
#   }
# }
validators <- list(
  cor_level=validate_numeric(min=0, max=1, min_warning=list(0.6, "correlation below 0.6 is departing from the spirit of this competition - to group things that likely contain redundant information")),
  k_splits=validate_numeric(min=2, max_warning=list(6, "these are a lot of splits...using this many splits with small data is probably unwise")),
  prevalence=validate_numeric(min=0, max=1),
  abundance=validate_numeric(min=0),
  lowest_level=validate_numeric(min=1, min_warning=list(2, "values below 2 may include an artificial taxonomic root")),
  max_level=validate_numeric(min=1, max=1000, max_warning=list(16, "you have many hierarchical levels, which may increase run time")),
  ncores=validate_numeric(min=1, max=parallel::detectCores()),
  nperm=validate_numeric(min=1, max=99999, max_warning=list(200, "this nperm value is high and will likely increase run time")),
  train_split=validate_numeric(min=0, max=1, min_warning=list(0.5, "a train test split below 50-50 is very unusual")),
  folds=validate_numeric(min=2, max_warning=list(11, "a value above 10 may result in very small splits")),
  tune_time=validate_numeric(min=0.1, max_warning=list(20, "spending excessive time tuning hyperparameters my not result in substaintal increases in accuracy")),
  permute=validate_numeric(min=1, max_warning=list(11, "you are about to permute the ML assessment pipeline more than 10 times, which is likely unnecessary")),
  seed=validate_numeric(min = -1 * .Machine$integer.max, max = .Machine$integer.max)
)

# Function to initialize parser for a program
# this takes in the data to make the parser but does not run it
initialize_parser <- function(version, program_name, description, argument_groups) {
  parser <- argparse::ArgumentParser(
    description=description,
    usage=paste(program_name, "[options] METADATA DATA"),
    formatter_class="type('CustomFormatter', (argparse.ArgumentDefaultsHelpFormatter, argparse.MetavarTypeHelpFormatter, argparse.RawTextHelpFormatter), {})"
  )

  # Common arguments for all programs
  parser$add_argument("METADATA", metavar="METADATA", type="character", help="path to metadata input (txt | tsv | csv)")
  parser$add_argument("DATA", metavar="DATA", type="character", help="path to input file from hierarchical data (i.e. hData data) (txt | tsv | csv)")

  parser$add_argument("-o", "--output_dir", type="character", metavar="<string>", default="outputs", help="Directory for the output files to be written. Defaults to a directory called 'outputs'")
  parser$add_argument("-v", "--version", action="version", version=version)
  parser$add_argument("--data_dir", type="character", metavar="<string>", default=".", help="Directory for MEATDATA, DATA, and output_dir, ignored if using absolute paths. Defaults to the current directory")
  parser$add_argument("--seed", type="numeric", metavar="<numeric>", default=default_seed(), help="Set the seed, if no value is provided, uses a random number from the range (-1 * 2^31, 2^31 - 1)")

  # add the arguments from the passed in argument_groups
  for (arg_group in argument_groups) {
    parser_group <- parser$add_argument_group(arg_group$name, arg_group$desc)
    for (arg in arg_group$args) {
      do.call(parser_group$add_argument, arg)
    }
  }

  return(parser)
}

# validate options loops through the validators and runs their validation functions
# if an error is returned by one of those functions, it will be cat-ed out and then quit is called
validate_options <- function(opts) {
  for (flag_to_validate in names(validators)) {
    # skip any validators not matching a flag
    if (is.null(opts[[flag_to_validate]])) next

    tryCatch({
      validators[[flag_to_validate]](flag_to_validate, opts[[flag_to_validate]], opts)
    }, error = function(e) {
      message(e$message)
      quit(status = 1)
    })
  }
}

# load the args for a program
# - loads version from env
# - initializes the parser with the desired argument groups
# - runs the parser and returns the arguments
load_args <- function(program_name, description, argument_groups) {
  # load version from the environment, defaulting to 0
  version <- Sys.getenv("TAXA_HFE_VERSION")
  if (version == "") {
    version <- "0"
  }

  # load the parser, including any arguments groups
  parser <- initialize_parser(version, program_name, description, argument_groups)

  # Parse the command-line arguments and return them
  opts <- parser$parse_args(commandArgs(TRUE))

  # run the validators against the parsed flags
  validate_options(opts)

  # extra handling for seed and files
  # because the validator has been run, it can be assumed that all options are safe for use

  # also normalize all input paths to the data_dir
  # will ignore the data_dir if the path links to valid file based on where the script is being run
  for (f in list("METADATA", "DATA", "output_dir")) {
    # lazyish check for abs path
    # change the file path to "data_dir / path" if the path doesn't start with "/"
    if (substr(opts[[f]], 0, 1) != "/") {
      opts[[f]] <- file.path(opts$data_dir, opts[[f]])
    }
  }

  # ensure the output directory does not have a trailing slash, and create it if it doesn't exist
  opts$output_dir <- gsub("/$", "", x = opts$output_dir)
  dir.create(opts$output_dir, showWarnings = FALSE)

  # set the seed from the flags
  set.seed(opts$seed)

  return(opts)
}

load_taxa_hfe_args <- function() {
  arg_groups <- list(argument_groups$taxa_hfe_base_args)

  return(load_args("taxa_hfe", "Hierarchical feature engineering (HFE) for feature reduction", arg_groups))
}

load_taxa_hfe_ml_args <- function() {
  arg_groups <- list(
    argument_groups$taxa_hfe_base_args,
    argument_groups$taxa_hfe_ml_args
  )

  return(load_args("taxa_hfe_ml", "Hierarchical feature engineering (HFE) with ML", arg_groups))
}
