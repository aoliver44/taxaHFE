#!/usr/bin/env Rscript

## SCRIPT: dietML.R ===================================================
## AUTHOR: Andrew Oliver
## DATE:   Nov, 1 2022
##
## PURPOSE: Run classification or regression ML
## the dietML.R script

## docker info =================================================================
#docker run --rm -v `pwd`:/home/docker -w /home/docker aoliver44/leakage_free_taxaHFE:latest

## suppress warnings
options(warn=-1)

## set working dir to /home for the docker container
setwd("/home/docker")

## load libraries ==============================================================
library(readr, quietly = T, verbose = F, warn.conflicts = F)
library(dplyr, quietly = T, verbose = F, warn.conflicts = F)
library(ranger, quietly = T, verbose = F, warn.conflicts = F)
suppressPackageStartupMessages(library(doParallel, quietly = T, verbose = F, warn.conflicts = F))
suppressPackageStartupMessages(library(parallel, quietly = T, verbose = F, warn.conflicts = F))
library(mikropml, quietly = T, verbose = F, warn.conflicts = F)
library(ggplot2, quietly = T, verbose = F, warn.conflicts = F)
library(fastshap, quietly = T, verbose = F, warn.conflicts = F)

## helper functions ============================================================
## Negate function ("not in"):
`%!in%` <- Negate(`%in%`)

## unregister hung-up parallel jobs
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

## list of supported models
models <- c("rf", "lasso", "ridge", "enet")

## check for inputs ============================================================

## quietly check to make sure /scripts wasn't overwritten
if (file.exists("/scripts/models/dietML_ranger_tidy.R") == FALSE) {
  stop("It appears you bind mounted docker to a virtual directory named /scripts. We
       need to use that folder. Please restart the docker image and use a different
       virtual directory name.")
}

## check for outdir and make if not there
if (dir.exists(paste0(dirname(opt$OUTPUT), "/ml_analysis")) == TRUE) {
  setwd(paste0(dirname(opt$OUTPUT), "/ml_analysis"))
} else {
  dir.create(path = paste0(dirname(opt$OUTPUT), "/ml_analysis"))
  setwd(paste0(dirname(opt$OUTPUT), "/ml_analysis"))
}

## read in inputs
train_data <- train_data_for_dietML
test_data <- test_data_for_dietML

## check for label
if ("feature_of_interest" %in% colnames(train_data) == FALSE & "feature_of_interest" %in% colnames(test_data) == FALSE) {
  stop(paste0("label not found in training AND testing data"))
} 

## check if classification was mis-specified
if (opt$feature_type == "factor") {
  type <<- "classification"
  if(length(levels(as.factor(metadata$feature_of_interest))) > 9)
  stop("You are trying to predict 10 or more classes. That is a bit much. Did you mean to do regression?")
} else {
  type <<- "regression"
}

## output all the parameters used ==============================================

cat("\n#########################\n")
cat("         DietML", "\n")
cat("#########################\n\n")

## run null (dummy) model ======================================================

cat("\n#########################\n")
cat("Running null model...", "\n")
cat("#########################\n\n")

source("/scripts/models/dietML_null_tidy.R")

## run chosen model ============================================================

## check if user input model
if (opt$model %!in% models) {
  cat("\n#########################\n")
  cat("ERROR: model not found", "\n")
  cat("Please choose one of the following models for --model ", "\n")
  print(as.data.frame(models))
  cat("#########################\n\n")
}

cat("\n#########################\n")
cat("Running model...", "\n")
cat("#########################\n\n")

## random forest
if (opt$model %in% c("ranger", "rf", "randomforest")) {
  source("/scripts/models/dietML_ranger_tidy.R")
}

## lasso/ridge models
if (opt$model %in% c("lasso", "ridge")) {
    source("/scripts/models/dietML_glmnet_tidy_ridge_lasso.R")
} 


## elastic net models
if (opt$model %in% c("enet", "elasticnet")) {
    source("/scripts/models/dietML_glmnet_tidy_enet.R")
} 

## VIP Plots ===================================================================
## For all:
# vip <- caret::varImp(object = training_fit)
# pdf(file = paste0(opt$outdir, "vip_plot.pdf"), width=15, height=5)
# plot(vip, top = pmin(NROW(vip$importance), 20))
# suppressMessages(dev.off())

## SHAP explanation ============================================================

if (opt$shap == TRUE) {
  
  cat("\n#########################\n")
  cat("Calculating Feature Importance", "\n")
  cat("#########################\n\n")
  
  if (opt$model %in% c("ranger", "rf", "randomforest")) {
    source("/scripts/utilities/shap_figures.R")
  } else if (opt$model %in% c("enet", "elasticnet", "lasso", "ridge")) {
    cat("You are not using a RF model. We will attempt to\n calculate VIP values instead.")
    source("/scripts/utilities/vip_basic.R")
  } else {
    cat("We are unable to calculate feature importances for\n the chosen model at this time.")
  }
}


## Done ========================================================================

save.image(file = paste0(dirname(opt$OUTPUT), "/ml_analysis/", "ML_r_workspace.rds"))

cat("\n#########################\n")
cat("Done! Results written to outdir.", "\n")
cat("#########################\n\n")
