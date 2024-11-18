#!/usr/bin/env Rscript
## v0.3.0a.8

## SCRIPT: vip_basic.R =========================================================
## AUTHOR: Andrew Oliver
## DATE:   Oct, 31 2023
##
## PURPOSE: fastshap and viz for tidymodels

## load libraries
library(vip, quietly = T, verbose = F, warn.conflicts = F)
library(ggplot2, quietly = T, verbose = F, warn.conflicts = F)
library(forcats, quietly = T, verbose = F, warn.conflicts = F)
library(tidymodels, quietly = T, verbose = F, warn.conflicts = F)

## catch the error
vip.error.occured <- FALSE

tryCatch( { if ((type == "classification") && (length(levels(as.factor(split_from_data_frame$data$feature_of_interest))) == 2) || (type == "regression")) {
  if (type == "regression") {
    gg_label <- opts$label
  } else {
    gg_label <- levels(as.factor(split_from_data_frame$data$feature_of_interest))[1]
  }
## Run on factor models: Full Data
final_fit <- best_tidy_workflow %>%
  tune::finalize_workflow(best_mod) %>%
  parsnip::fit(data = split_from_data_frame$data) %>%
  workflowsets::extract_fit_parsnip()

importance_plot_full <- final_fit %>%
  vip::vi(lambda = best_mod$penalty) %>%
  dplyr::mutate(Variable = fct_reorder(Variable, Importance)) %>%
  dplyr::slice_head(n = 25) %>%
  ggplot2::ggplot(aes(x = Importance, y = Variable, fill = Sign)) +
  geom_col() +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = NULL) +
  ggtitle(label = paste0("Predicting: ", gg_label), subtitle = "using all data") + 
  theme_minimal()

ggplot2::ggsave(plot = importance_plot_full, 
                filename = paste0(dirname(opts$OUTPUT), "/ml_analysis/", "importance_plot_full.pdf"),
                width = pmax((0.1 * max(nchar(colnames(input)))), 6), height = 4.5, units = "in")

## Run on factor models: Train Data
final_fit <- best_tidy_workflow %>%
  tune::finalize_workflow(best_mod) %>%
  parsnip::fit(data = train) %>%
  workflowsets::extract_fit_parsnip()

importance_plot_train <- final_fit %>%
  vip::vi(lambda = best_mod$penalty) %>%
  dplyr::mutate(Variable = fct_reorder(Variable, Importance)) %>%
  dplyr::slice_head(n = 25) %>%
  ggplot2::ggplot(aes(x = Importance, y = Variable, fill = Sign)) +
  geom_col() +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = NULL) +
  ggtitle(label = paste0("Predicting: ", gg_label), subtitle = "using train data") + 
  theme_minimal()

ggplot2::ggsave(plot = importance_plot_train, 
                filename = paste0(dirname(opts$OUTPUT), "/ml_analysis/", "importance_plot_train.pdf"),
                width = pmax((0.1 * max(nchar(colnames(train)))), 6), height = 4.5, units = "in")

## Run on factor models: Test Data
final_fit <- best_tidy_workflow %>%
  tune::finalize_workflow(best_mod) %>%
  parsnip::fit(data = test) %>%
  workflowsets::extract_fit_parsnip()

importance_plot_test <- final_fit %>%
  vip::vi(lambda = best_mod$penalty) %>%
  dplyr::mutate(Variable = fct_reorder(Variable, Importance)) %>%
  dplyr::slice_head(n = 25) %>%
  ggplot2::ggplot(aes(x = Importance, y = Variable, fill = Sign)) +
  geom_col() +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = NULL) +
  ggtitle(label = paste0("Predicting: ", gg_label), subtitle = "using test data") + 
  theme_minimal()

ggplot2::ggsave(plot = importance_plot_test, 
                filename = paste0(dirname(opts$OUTPUT), "/ml_analysis/", "importance_plot_test.pdf"),
                width = pmax((0.1 * max(nchar(colnames(test)))), 6), height = 4.5, units = "in")
}
  
}, error = function(e) {vip.error.occured <<- TRUE} )

if (vip.error.occured == TRUE) {
  cat("\n#########################\n")
  cat("ERROR: Could not complete feature importance anlaysis.", "\n")
  cat("#########################\n\n")
}
  