#!/usr/bin/env Rscript
## v0.3.0a.8

## SCRIPT: diet_ml_null_tidy.R ===================================================
## AUTHOR: Andrew Oliver
## DATE:   Mar, 20 2023
##
## PURPOSE: NULL model for tidymodels

## helper functions and vars ===================================================

## suppress warnings
options(warn=-1)

## Negate function ("not in"):
`%!in%` <- Negate(`%in%`)

## set seed
set.seed(as.numeric(opts$seed))

## load libraries ==============================================================

library(tidyr, quietly = T, verbose = F, warn.conflicts = F)
library(ggplot2, quietly = T, verbose = F, warn.conflicts = F)
suppressPackageStartupMessages(library(tidymodels, quietly = T, verbose = F, warn.conflicts = F))

## create results df ===========================================================

if (opts$feature_type == "factor") {
  results_df <- data.frame(seed = numeric(), bal_accuracy = numeric(), f_meas = numeric(), accuracy = numeric(), stringsAsFactors = F)
} else if (opts$feature_type == "numeric") {
  results_df <- data.frame(seed = numeric(), mae = numeric(), rmse = numeric(), ccc = numeric(), stringsAsFactors = F)
}

## interate over null model ====================================================

  
if (opts$feature_type == "factor") {
  df_loop_results <- data.frame(truth = character(), estimate = character(), stringsAsFactors = F)
} else if (opts$feature_type == "numeric") {
  df_loop_results <- data.frame(truth = numeric(), estimate = numeric(), stringsAsFactors = F)
}

## set initial test-train split
train <- train_data
test  <- test_data

## recipe ======================================================================

## specify recipe (this is like the pre-process work)
diet_ml_recipe <- 
  recipes::recipe(feature_of_interest ~ ., data = train) %>% 
  recipes::update_role(tidyr::any_of("subject_id"), new_role = "ID") %>%
  recipes::step_dummy(recipes::all_nominal_predictors()) %>%
  recipes::step_corr(all_numeric_predictors(), threshold = as.numeric(opts$cor_level)) %>%
  recipes::step_zv(all_predictors())


## ML engine ===================================================================

## specify ML model and engine 
initial_mod <- null_model() %>% 
  set_engine("parsnip") %>% 
  set_mode(type) %>% 
  translate()

## workflow ====================================================================

## define workflow
diet_ml_wflow <- 
  workflows::workflow() %>% 
  workflows::add_model(initial_mod) %>% 
  workflows::add_recipe(diet_ml_recipe)  


## fit model ==============================================================

## fit to test data
final_res <- parsnip::fit(diet_ml_wflow, test)

df_loop_results <- add_row(df_loop_results, truth = test$feature_of_interest)
df_loop_results$estimate <- final_res$fit$fit$fit$value

if (opts$feature_type == "factor") {
  df_loop_results$estimate <- factor(x = df_loop_results$estimate, levels = levels(as.factor(df_loop_results$truth)))
  results_df <- results_df %>% 
    tibble::add_row(., bal_accuracy = 
                      yardstick::bal_accuracy_vec(truth = as.factor(df_loop_results$truth), 
                                     estimate = as.factor(df_loop_results$estimate), 
                                     data = df_loop_results), 
                    accuracy = 
                      yardstick::accuracy_vec(truth = as.factor(df_loop_results$truth), 
                                              estimate = as.factor(df_loop_results$estimate), 
                                              data = df_loop_results), 
                    f_meas = 
                      yardstick::f_meas_vec(truth = as.factor(df_loop_results$truth), 
                                              estimate = as.factor(df_loop_results$estimate), 
                                              data = df_loop_results),
                    seed = opts$seed)
} else if (opts$feature_type == "numeric") {
  results_df <- results_df %>% 
    tibble::add_row(., mae = 
                      yardstick::mae_vec(truth = df_loop_results$truth, 
                                     estimate = df_loop_results$estimate, 
                                     data = df_loop_results), 
                    rmse = 
                      yardstick::rmse_vec(truth = df_loop_results$truth, 
                                     estimate = df_loop_results$estimate, 
                                     data = df_loop_results),
                    ccc = yardstick::ccc_vec(truth = df_loop_results$truth, 
                                          estimate = df_loop_results$estimate, 
                                          data = df_loop_results),
                    seed = opts$seed)
  
}
  


## write df ====================================================================

## write table of results to file
readr::write_csv(x = results_df, file ="dummy_model_results.csv", 
                 append = T, col_names = !file.exists("dummy_model_results.csv"))

## show the final results
cat("Performance of NULL model:", "\n")
cat("Label: ", opts$label, "\n")
print(results_df %>% dplyr::select(-seed) %>% dplyr::summarise_all(., ~mean(.x)))

