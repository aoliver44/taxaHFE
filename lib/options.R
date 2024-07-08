## holds the parser for the flags of taxaHFE
suppressPackageStartupMessages(library(argparse, quietly = T, verbose = F, warn.conflicts = F))

# in order to reduce code and specify different descriptions, the arg_parser()
# function requires you to specify the program_list_option as a numeric.
# program_list_option:
# 1. taxaHFE
# 2. taxaHFE-ML

arg_parser <- function(version, program_list_option) {
  parser <- argparse::ArgumentParser(
    description=switch(program_list_option,
      {'Hierarchical feature engineering (HFE) for the reduction of features with respects to a factor or regressor'},
      {'Hierarchical feature engineering (HFE) for the reduction of features with respects to a factor or regressor, using a train-test split and machine learning to identify features of hierarchical features of interest'} 
      ),
    usage=paste(switch(program_list_option, {'taxaHFE'},{'taxaHFE-ML'}), "[options] METADATA DATA OUTPUT"),
    ## this is a grouping of built in argparse formatter to make the help printout be more clear
    ## ArgumentDefaultsHelpFormatter - displays the default next to the help text
    ## MetavarTypeHelpFormatter - displays the expected type for the flag next to the flag names
    ## RawTextHelpFormatter - less line-wrapping in the printout
    ## https://docs.python.org/3.9/library/argparse.html (library uses OS python3 version)
    formatter_class="type('CustomFormatter', (argparse.ArgumentDefaultsHelpFormatter, argparse.MetavarTypeHelpFormatter, argparse.RawTextHelpFormatter), {})"
  )
  
  parser$add_argument('METADATA', metavar='METADATA', type="character", help="path to metadata input (txt | tsv | csv)")
  parser$add_argument('DATA', metavar='DATA', type="character", help="path to input file from hierarchical data (i.e. hData data) (txt | tsv | csv)")
  parser$add_argument('OUTPUT', metavar='OUTPUT', type="character", help="output file name (csv)")
  
  taxahfe_base_args=parser$add_argument_group('TaxaHFE arguments', 'Options to pass to TaxaHFE')
  taxahfe_base_args$add_argument('-s', '--subject_identifier', type="character", metavar="<string>", default='subject_id', help='Metadata column name containing subject IDs')
  taxahfe_base_args$add_argument('-l', '--label', type="character", metavar="<string>", default='feature_of_interest', help='Metadata column name of interest for ML')
  taxahfe_base_args$add_argument('-t', '--feature_type', type="character", metavar="<string>", default='factor', help='Is the ML label a factor or numeric')
  taxahfe_base_args$add_argument('-a', '--abundance', type="numeric", metavar="<numeric>", default='0', help='Minimum mean abundance of feature')
  taxahfe_base_args$add_argument('-p', '--prevalence', type="numeric", metavar="<numeric>", default='0.01', help='Minimum prevalence of feature')
  taxahfe_base_args$add_argument('-L', '--lowest_level', type="integer", metavar="<numeric>", default='2', help='Most general level allowed to compete')
  taxahfe_base_args$add_argument('-m', '--max_level', type="integer", metavar="<numeric>", default='1000', help='How many hierarchical levels should be allowed to compete')
  taxahfe_base_args$add_argument('-c', '--cor_level', type="numeric", metavar="<numeric>", default='0.95', help='Initial pearson correlation filter')
  taxahfe_base_args$add_argument('-d', '--disable_super_filter', action="store_true", help='Disable running of the super filter (final forest competition)')
  taxahfe_base_args$add_argument('-w', '--write_old_files', action="store_true", help='Write individual level files and old HFE files')
  taxahfe_base_args$add_argument('-W', '--write_flattened_tree', action="store_true", help='Write a compressed backup of the entire competed tree')
  taxahfe_base_args$add_argument('-D', '--write_both_outputs', action="store_true", help='Write an output for pre and post super filter results, overridden by --disable_super_filter')
  taxahfe_base_args$add_argument('--nperm', type="integer", metavar="<numeric>", default='40', help='Number of RF permutations')
  parser$add_argument('-n', '--ncores', type="integer", metavar="<numeric>", default='2', help='Number of cpu cores to use')
  parser$add_argument('--seed', type="numeric", metavar="<numeric>", help='Set a random numeric seed. If None, defaults to use system time')
  parser$add_argument('-v', '--version', action='version', version=version)
  return(parser)
}

taxaHFE_ML_arg_parser <- function(version, program_list_option) {
  parser <- arg_parser(version, program_list_option)
  
  dietml_base_args=parser$add_argument_group('DietML arguments', 'Options to pass to DietML for machine learning and SHAP analysis of TaxaHFE features')
  dietml_base_args$add_argument('--train_split', type="numeric", metavar="<numeric>", default='0.8', help='Percentage of samples should be used in training')
  dietml_base_args$add_argument('--model', type="character", metavar="<string>", default='rf', choices=c("rf","enet", "none"), help='Percentage of samples should be used in training')
  dietml_base_args$add_argument('--folds', type="numeric", metavar="<numeric>", default='10', help='Number of CV folds to tune with')
  dietml_base_args$add_argument('--metric', type="character", metavar="<string>", default='bal_accuracy', choices=c("roc_auc", "bal_accuracy", "accuracy", "mae", "rmse", "rsq", "kap", 
                                                                                                                    "f_meas", "ccc"), help='Metric would you like to optimize in training')
  dietml_base_args$add_argument('--tune_length', type="numeric", metavar="<numeric>", default='80', help='Number of hyperparameter combinations to sample')
  dietml_base_args$add_argument('--tune_time', type="numeric", metavar="<numeric>", default='10', help='Length of time hyperparameter search runs')
  dietml_base_args$add_argument('--tune_stop', type="numeric", metavar="<numeric>", default='10', help='Number of HP interations to let pass without a metric improvement')
  dietml_base_args$add_argument('--shap', action="store_true", help='Attempt to calcualte shap values?')
  return(parser)
}

load_args <- function(version, program_list_option) {
  parser <- switch(program_list_option, 
                   {arg_parser(version, program_list_option)},
                   {taxaHFE_ML_arg_parser(version, program_list_option)}
                   )
  return(parser$parse_args(commandArgs(TRUE)))
}

