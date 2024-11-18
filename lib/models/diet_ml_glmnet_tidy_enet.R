#!/usr/bin/env Rscript
## v0.3.0a.8

## SCRIPT: diet_ml_glmnet_tidy_enet.R ===================================================
## AUTHOR: Andrew Oliver
## DATE:   Jan, 30 2023
##
## PURPOSE: Glmnet model for tidymodels

## helper functions and vars ===================================================

## suppress warnings
options(warn=-1)

## Negate function ("not in"):
`%!in%` <- Negate(`%in%`)

## unregister hung-up parallel jobs
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

## set seed
set.seed(as.numeric(opts$seed))

## load libraries ==============================================================

library(mikropml, quietly = T, verbose = F, warn.conflicts = F)
library(tidyr, quietly = T, verbose = F, warn.conflicts = F)
library(ggplot2, quietly = T, verbose = F, warn.conflicts = F)
suppressPackageStartupMessages(library(tidymodels, quietly = T, verbose = F, warn.conflicts = F))
suppressPackageStartupMessages(library(glmnet, quietly = T, verbose = F, warn.conflicts = F))


## resample strategy ===========================================================

## set initial test-train split
train <- train_data
test  <- test_data
split_from_data_frame <- make_splits(
  x = train,
  assessment = test
)

## set resampling scheme
folds <- rsample::vfold_cv(train, v = as.numeric(opts$folds), strata = label, repeats = 3)

## recipe ======================================================================

## specify recipe (this is like the pre-process work)
diet_ml_recipe <- 
  recipes::recipe(feature_of_interest ~ ., data = train) %>% 
  recipes::update_role(tidyr::any_of("subject_id"), new_role = "ID") %>% 
  recipes::step_dummy(recipes::all_nominal_predictors()) %>%
  recipes::step_corr(recipes::all_numeric_predictors(), threshold = as.numeric(opts$cor_level), use = "everything") %>%
  recipes::step_zv(recipes::all_predictors())

## ML engine ===================================================================

## specify ML model and engine 
if (type == "classification") {
  initial_mod <- parsnip::logistic_reg(mode = "classification", 
                                       penalty = tune(),
                                       mixture = tune()) %>%
    parsnip::set_engine("glmnet")
} else {
  initial_mod <- parsnip::linear_reg(mode = "regression", 
                                       penalty = tune(),
                                       mixture = tune()) %>%
    parsnip::set_engine("glmnet")
}

initial_mod %>% parsnip::translate()

## workflow ====================================================================

## define workflow
diet_ml_wflow <- 
  workflows::workflow() %>% 
  workflows::add_model(initial_mod) %>% 
  workflows::add_recipe(diet_ml_recipe)  

## set up parallel jobs ========================================================
## remove any doParallel job setups that may have
## unneccessarily hung around
unregister_dopar()

## register parallel cluster
cl <- parallel::makePSOCKcluster(as.numeric(opts$ncores))
doParallel::registerDoParallel(cl)

## hyperparameters =============================================================

## define the hyper parameter set
diet_ml_param_set <- parsnip::extract_parameter_set_dials(diet_ml_wflow)

## set up hyper parameter search
if (type == "classification") {
  search_res <-
    diet_ml_wflow %>% 
    tune::tune_bayes(
      resamples = folds,
      # To use non-default parameter ranges
      param_info = diet_ml_param_set,
      # Generate five at semi-random to start
      initial = 5,
      iter = opts$tune_length,
      # How to measure performance?
      metrics = yardstick::metric_set(bal_accuracy, roc_auc, accuracy, kap),
      control = tune::control_bayes(no_improve = as.numeric(opts$tune_stop),
                                    uncertain = 5,
                                    verbose = FALSE,
                                    parallel_over = "resamples",
                                    time_limit = as.numeric(opts$tune_time),
                                    seed = as.numeric(opts$seed))
    )
  
} else if (type == "regression") {
  search_res <-
    diet_ml_wflow %>% 
    tune::tune_bayes(
      resamples = folds,
      # To use non-default parameter ranges
      param_info = diet_ml_param_set,
      # Generate five at semi-random to start
      initial = 5,
      iter = opts$tune_length,
      # How to measure performance?
      metrics = yardstick::metric_set(mae, rmse, rsq, ccc),
      control = tune::control_bayes(no_improve = as.numeric(opts$tune_stop),
                                    uncertain = 5,
                                    verbose = FALSE,
                                    parallel_over = "resamples",
                                    time_limit = as.numeric(opts$tune_time),
                                    seed = as.numeric(opts$seed))
    )
}

search_res %>% tune::show_best(opts$metric)

## stop parallel jobs
parallel::stopCluster(cl)
## remove any doParallel job setups that may have
## unneccessarily hung around
unregister_dopar()

## fit best model ==============================================================

## get the best parameters from tuning
best_mod <- 
  search_res %>% 
  tune::select_best(metric = opts$metric)

## create the last model based on best parameters
if (type == "classification") {
  last_best_mod <- 
    parsnip::logistic_reg(mode = "classification", penalty = best_mod$penalty, mixture = best_mod$mixture) %>% 
    parsnip::set_engine("glmnet") %>% 
    parsnip::set_mode(type)
} else {
  last_best_mod <- 
    parsnip::linear_reg(mode = "regression", penalty = best_mod$penalty, mixture = best_mod$mixture) %>% 
    parsnip::set_engine("glmnet") %>% 
    parsnip::set_mode(type)
}

## update workflow with best model
best_tidy_workflow <- 
  diet_ml_wflow %>% 
  workflows::update_model(last_best_mod)

## fit to test data
if (type == "classification") {
  final_res <- tune::last_fit(best_tidy_workflow, split_from_data_frame, 
                              metrics = yardstick::metric_set(bal_accuracy, 
                                                              roc_auc, accuracy, 
                                                              kap, f_meas))
} else if (type == "regression") {
  final_res <- tune::last_fit(best_tidy_workflow, split_from_data_frame, 
                              metrics = yardstick::metric_set(mae, rmse, rsq, 
                                                              ccc))
}


cat("\n################\n")
cat("RESULTS:", "\n")
cat("##################\n\n")

## show the final results
cat("Performance of test set:", "\n")
cat("File: ", opts$input, "\n")
cat("Label: ", opts$label, "\n")
cat("Model: ", opts$model, "\n")
print(workflowsets::collect_metrics(final_res))

## merge null results with trained results and write table
null_results <- results_df %>% 
  dplyr::select(., -seed) %>% 
  summarise_all(., mean) %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = ".metric") %>% 
  dplyr::rename(., "null_model_avg" = 2)
full_results <- merge(workflowsets::collect_metrics(final_res), null_results, by = ".metric", all = T)
full_results$seed <- opts$seed

## write final results to file or append if file exists
readr::write_csv(x = full_results, file = paste0(opts$outdir, "ml_results.csv"), append = T, col_names = !file.exists(paste0(opts$outdir, "ml_results.csv")))


## graphs ======================================================================

hyperpar_perf_plot <- autoplot(search_res, type = "performance")
ggplot2::ggsave(plot = hyperpar_perf_plot, filename = paste0(opts$outdir, "training_performance.pdf"), width = 7, height = 2.5, units = "in")

hyperpar_tested_plot <- autoplot(search_res, type = "parameters") + 
  labs(x = "Iterations", y = NULL)
ggplot2::ggsave(plot = hyperpar_tested_plot, filename = paste0(opts$outdir, "hyperpars_tested.pdf"), width = 7, height = 2.5, units = "in")

