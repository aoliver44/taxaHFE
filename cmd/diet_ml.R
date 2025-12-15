source("lib/options.R")

# load flags
# to use this code line-by-line in the Rstudio context, commandArgs can be overloaded to specify the desired flags
# ex. commandArgs <- function(x) { c("example_inputs/metadata.txt", "example_inputs/microbiome_data.txt", "-o", "example_outputs", "-s", "Sample", "-l", "Category", "-L", "3", "-n", "2", "--seed", "42", "--train_split", "0.8", "--tune_time", "0", "-W", "--parallel_workers", "2", "--shap") }
# these will be used by the argparser
opts <- load_diet_ml_args()

# TODO: the rest
print(opts)