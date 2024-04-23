#!/usr/bin/env Rscript
## v0.3.0a.8

## SCRIPT: shap_figures.R ======================================================
## AUTHOR: Andrew Oliver
## DATE:   Jan, 30 2023
##
## PURPOSE: fastshap and viz for tidymodels

## load libraries

library(fastshap, quietly = T, verbose = F, warn.conflicts = F)
library(shapviz, quietly = T, verbose = F, warn.conflicts = F)
library(ggplot2, quietly = T, verbose = F, warn.conflicts = F)
library(tidymodels, quietly = T, verbose = F, warn.conflicts = F)
library(recipes, quietly = T, verbose = F, warn.conflicts = F)

shap.error.occured <- FALSE

tryCatch( { if (length(levels(as.factor(split_from_data_frame$data$feature_of_interest))) == 2) {
  
  #####################
  ## FULL DATA 1
  #####################
  
  ## Prediction wrapper
  pfun <- function(object, newdata) {
    predict(object, data = newdata)$predictions[, levels(as.factor(split_from_data_frame$data$feature_of_interest))[1]]
  }
  
  ## pull model out of workflow
  best_workflow <- best_tidy_workflow %>%
    parsnip::fit(split_from_data_frame$data)
  best_workflow_mod <- workflows::extract_fit_parsnip(best_workflow)
  
  ## pull out data
  shap_data_full <- recipes::prep(dietML_recipe, split_from_data_frame$data) %>% 
    recipes::juice() %>% 
    dplyr::select(-feature_of_interest, -dplyr::any_of(opt$subject_identifier)) %>% 
    as.matrix()
  
  ## explain with fastshap
  shap_explainations_full <- fastshap::explain(best_workflow_mod$fit, X = shap_data_full, pred_wrapper = pfun, nsim = 100, adjust = TRUE)
  
  ## make shap viz object
  sv_full <- shapviz::shapviz(shap_explainations_full, X = shap_data_full)
  
  ## make shap plot
  importance_plot_full_1 <- shapviz::sv_importance(sv_full, kind = "bee", show_numbers = TRUE, bee_width = 0.2, max_display = 10) + 
    ggtitle(label = paste0("SHAP: ", opt$label, " (full data)")) + 
    labs(x = paste0("predictive of ", levels(as.factor(input$feature_of_interest))[2], " < SHAP > ", "predictive of ", levels(as.factor(input$feature_of_interest))[1])) + 
    theme_bw(base_size = 14)
  ggplot2::ggsave(plot = importance_plot_full_1, filename = paste0(dirname(opt$OUTPUT), "/ml_analysis/", "importance_plot_full_1.pdf"), width = pmax((0.1 * max(nchar(colnames(sv_full$X)))), 6), height = 4.5, units = "in")
  
  #####################
  ## TRAIN DATA 1
  #####################
  
  ## Prediction wrapper: first level (ie if levels are high, low, this is high)
  pfun <- function(object, newdata) {
    predict(object, data = newdata)$predictions[, levels(as.factor(split_from_data_frame$data$feature_of_interest))[1]]
  }
  
  ## pull model out of workflow
  best_workflow <- best_tidy_workflow %>%
    parsnip::fit(train)
  best_workflow_mod <- workflows::extract_fit_parsnip(best_workflow)
  
  ## pull out data
  shap_data_train <- recipes::prep(dietML_recipe, train) %>% 
    recipes::juice() %>% 
    dplyr::select(-feature_of_interest, -dplyr::any_of(opt$subject_identifier)) %>% 
    as.matrix()
  
  ## explain with fastshap
  shap_explainations_train <- fastshap::explain(best_workflow_mod$fit, X = shap_data_train, pred_wrapper = pfun, nsim = 100, adjust = TRUE)

  ## make shap viz object
  sv_train <- shapviz::shapviz(shap_explainations_train, X = shap_data_train)
  
  ## make shap plot
  importance_plot_train_1 <- shapviz::sv_importance(sv_train, kind = "bee", show_numbers = TRUE, bee_width = 0.2, max_display = 10) + 
    ggtitle(label = paste0("SHAP: ", opt$label, " (train data)")) + 
    labs(x = paste0("predictive of ", levels(as.factor(input$feature_of_interest))[2], " < SHAP > ", "predictive of ", levels(as.factor(input$feature_of_interest))[1])) + 
    theme_bw(base_size = 14)
  ggplot2::ggsave(plot = importance_plot_train_1, filename = paste0(dirname(opt$OUTPUT), "/ml_analysis/", "importance_plot_train_1.pdf"), width = pmax((0.1 * max(nchar(colnames(sv_train$X)))), 6), height = 4.5, units = "in")
  
  #####################
  ## TEST DATA 1
  #####################
  
  ## Prediction wrapper: second level (ie if levels are high, low, this is low)
  pfun <- function(object, newdata) {
    predict(object, data = newdata)$predictions[, levels(as.factor(split_from_data_frame$data$feature_of_interest))[1]]
  }
  
  ## pull model out of workflow
  best_workflow <- best_tidy_workflow %>%
    parsnip::fit(test)
  best_workflow_mod <- workflows::extract_fit_parsnip(best_workflow)
  
  ## pull out data
  shap_data_test<- recipes::prep(dietML_recipe, test) %>% 
    recipes::juice() %>% 
    dplyr::select(-feature_of_interest, -dplyr::any_of(opt$subject_identifier)) %>% 
    as.matrix()
  
  ## explain with fastshap
  shap_explainations_test <- fastshap::explain(best_workflow_mod$fit, X = shap_data_test, pred_wrapper = pfun, nsim = 100, adjust = TRUE)
  
  ## make shap viz object
  sv_test <- shapviz::shapviz(shap_explainations_test, X = shap_data_test)
  
  ## make shap plot
  importance_plot_test_1 <- shapviz::sv_importance(sv_test, kind = "bee", show_numbers = TRUE, bee_width = 0.2, max_display = 10) + 
    ggtitle(label = paste0("SHAP: ", opt$label, " (test data)")) + 
    labs(x = paste0("predictive of ", levels(as.factor(input$feature_of_interest))[2], " < SHAP > ", "predictive of ", levels(as.factor(input$feature_of_interest))[1])) + 
    theme_bw(base_size = 14)
  ggplot2::ggsave(plot = importance_plot_test_1, filename = paste0(dirname(opt$OUTPUT), "/ml_analysis/", "importance_plot_test_1.pdf"), width = pmax((0.1 * max(nchar(colnames(sv_test$X)))), 6), height = 4.5, units = "in")
  
  } 
  
  if ((type == "regression") && (opt$model == "rf")) {
    
    ## Prediction wrapper
    pfun <- function(object, newdata) {
      predict(object, data = newdata)$predictions
    }
    
    ## pull model out of workflow
    best_workflow <- best_tidy_workflow %>%
      parsnip::fit(split_from_data_frame$data)
    best_workflow_mod <- workflows::extract_fit_parsnip(best_workflow)
    
    ## pull out data
    shap_data_full <- recipes::prep(dietML_recipe, split_from_data_frame$data) %>% 
      recipes::juice() %>% 
      dplyr::select(-feature_of_interest, -dplyr::any_of(opt$subject_identifier)) %>% 
      as.matrix()
    
    ## explain with fastshap
    shap_explainations_full <- fastshap::explain(best_workflow_mod$fit, X = shap_data_full, pred_wrapper = pfun, nsim = 100, adjust = TRUE)
    
    ## make shap viz object
    sv_full <- shapviz::shapviz(shap_explainations_full, X = shap_data_full)
    
    ## make shap plot
    importance_plot_full <- shapviz::sv_importance(sv_full, kind = "bee", show_numbers = TRUE, bee_width = 0.2, max_display = 10) + 
      ggtitle(label = paste0("SHAP: ", opt$label, " (full data)")) + theme_bw(base_size = 14)
    ggplot2::ggsave(plot = importance_plot_full, filename = paste0(dirname(opt$OUTPUT), "/ml_analysis/", "importance_plot_full.pdf"), width = pmax((0.1 * max(nchar(colnames(sv_full$X)))), 6), height = 4.5, units = "in")
    
    #######################
    
    ## Prediction wrapper
    pfun <- function(object, newdata) {
      predict(object, data = newdata)$predictions
    }
    
    ## pull model out of workflow
    best_workflow <- best_tidy_workflow %>%
      parsnip::fit(train)
    best_workflow_mod <- workflows::extract_fit_parsnip(best_workflow)
    
    ## pull out data
    shap_data_train <- recipes::prep(dietML_recipe, train) %>% 
      recipes::juice() %>% 
      dplyr::select(-feature_of_interest, -dplyr::any_of(opt$subject_identifier)) %>% 
      as.matrix()
    
    ## explain with fastshap
    shap_explainations_train <- fastshap::explain(best_workflow_mod$fit, X = shap_data_train, pred_wrapper = pfun, nsim = 100, adjust = TRUE)
    
    ## make shap viz object
    sv_train <- shapviz::shapviz(shap_explainations_train, X = shap_data_train)
    
    ## make shap plot
    importance_plot_train <- shapviz::sv_importance(sv_train, kind = "bee", show_numbers = TRUE, bee_width = 0.2, max_display = 10) + 
      ggtitle(label = paste0("SHAP: ", opt$label, " (train)")) + theme_bw(base_size = 14)
    ggplot2::ggsave(plot = importance_plot_train, filename = paste0(dirname(opt$OUTPUT), "/ml_analysis/", "importance_plot_train.pdf"), width = pmax((0.1 * max(nchar(colnames(sv_train$X)))), 6), height = 4.5, units = "in")
    
    ####################
    
    ## Prediction wrapper
    pfun <- function(object, newdata) {
      predict(object, data = newdata)$predictions
    }
    
    ## pull model out of workflow
    best_workflow <- best_tidy_workflow %>%
      parsnip::fit(test)
    best_workflow_mod <- workflows::extract_fit_parsnip(best_workflow)
    
    ## pull out data
    shap_data_test <- recipes::prep(dietML_recipe, test) %>% 
      recipes::juice() %>% 
      dplyr::select(-feature_of_interest, -dplyr::any_of(opt$subject_identifier)) %>% 
      as.matrix()
    
    ## explain with fastshap
    shap_explainations_test <- fastshap::explain(best_workflow_mod$fit, X = shap_data_test, pred_wrapper = pfun, nsim = 100, adjust = TRUE)
    
    ## make shap viz object
    sv_test <- shapviz::shapviz(shap_explainations_test, X = shap_data_test)
    
    ## make shap plot
    importance_plot_test <- shapviz::sv_importance(sv_test, kind = "bee", show_numbers = TRUE, bee_width = 0.2, max_display = 10) + 
      ggtitle(label = paste0("SHAP: ", opt$label, " (test)")) + theme_bw(base_size = 14)
    ggplot2::ggsave(plot = importance_plot_test, filename = paste0(dirname(opt$OUTPUT), "/ml_analysis/", "importance_plot_test.pdf"), width = pmax((0.1 * max(nchar(colnames(sv_test$X)))), 6), height = 4.5, units = "in")
    
    ####################
    
  }

}, error = function(e) {shap.error.occured <<- TRUE} )

if (shap.error.occured == TRUE) {
  cat("\n#########################\n")
  cat("ERROR: Could not complete SHAP anlaysis.", "\n")
  cat("#########################\n\n")
}
