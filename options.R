## holds the parser for the flags of taxaHFE
library(argparse, quietly = T, verbose = F, warn.conflicts = F)

arg_parser <- function(version) {
  parser <- argparse::ArgumentParser(
    description='Hierarchical feature engineering (HFE) for the reduction of features with respects to a factor or regressor',
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
  
  parser$add_argument('-v', '--version', action='version', version=version)
  parser$add_argument('-s', '--subject_identifier', type="character", default='subject_id', help='Metadata column name containing subject IDs')
  parser$add_argument('-l', '--label', type="character", default='cluster', help='Metadata column name of interest for ML')
  parser$add_argument('-t', '--feature_type', type="character", default='factor', help='Is the ML label a factor or numeric')
  parser$add_argument('-a', '--abundance', type="numeric", default='0.0001', help='Minimum mean abundance of feature')
  parser$add_argument('-p', '--prevalence', type="numeric", default='0.01', help='Minimum prevalence of feature')
  parser$add_argument('-L', '--lowest_level', type="integer", default='2', help='Most general level allowed to compete')
  parser$add_argument('-m', '--max_depth', type="integer", default='1000', help='How many hierarchical levels should be allowed to compete')
  parser$add_argument('-c', '--cor_level', type="numeric", default='0.95', help='Initial pearson correlation filter')
  parser$add_argument('-n', '--ncores', type="integer", default='2', help='Number of cpu cores to use')
  parser$add_argument('-d', '--disable_super_filter', action="store_true", help='Disable running of the super filter (final forest competition)')
  parser$add_argument('-w', '--write_old_files', action="store_true", help='Write individual level files and old HFE files')
  parser$add_argument('-W', '--write_flattened_tree', action="store_true", help='Write a compressed backup of the entire competed tree')
  parser$add_argument('-D', '--write_both_outputs', action="store_true", help='Write an output for pre and post super filter results, overridden by --disable_super_filter')
  parser$add_argument('--nperm', type="integer", default='40', help='Number of RF permutations')
  parser$add_argument('--seed', type="numeric", help='Set a random numeric seed. If None, defaults to use system time')
  
  return(parser)
}

load_args <- function(version) {
  parser <- arg_parser(version)
  return(parser$parse_args(commandArgs(TRUE)))
}
