## DIETML FUNCTIONS

## libraries  ==================================================================
source("lib/requirements.R")

## helper functions ============================================================

## Negate function ("not in"):
`%!in%` <- Negate(`%in%`)

## suppress warnings
options(warn = -1)

## run dietML based on diet_ml_input_df ========================================
run_diet_ml <- function(diet_ml_inputs, metadata, n_repeat, feature_type, seed, train, test, model, program, random_effects, folds, cv_repeats, ncores, parallel_workers, tune_length, tune_stop, tune_time, metric, label, output, shap) {
  input_df <- extract_attributes(diet_ml_inputs)
  
  ## create a number of random seeds, which are used across all programs 
  ## run, ie taxahfe, taxahfe-ml - they need to be run across the same seeds!
  ## to make sure these are created determistically based on set seed, make
  ## sure its set before (might not matter?)
  set.seed(seed)
  random_seeds <- sample(1:100000000, replace = F, size = n_repeat)
  ## make sure the first seed used is always the one specified. 
  random_seeds[1] <- seed
  
  for (seed in random_seeds) {
    for (dML_input in unique(input_df[["general_name"]])) {
      
      train_data <- diet_ml_inputs[[input_df %>% 
                                      dplyr::filter(., general_name == dML_input & train_test_attr == "train") %>% 
                                      dplyr::pull(name)]]
      
      test_data <- diet_ml_inputs[[input_df %>% 
                                     dplyr::filter(., general_name == dML_input & train_test_attr == "test") %>% 
                                     dplyr::pull(name)]]
      
      ## keep track of what method is being passed to dietML
      ## this gets printed in the results file
      program <- dML_input
      pass_to_dietML(train = train_data, 
                     test = test_data,
                     metadata = metadata,
                     program = program, 
                     model = model, 
                     seed = seed, 
                     random_effects = random_effects, 
                     folds = folds,
                     cv_repeats = cv_repeats,
                     ncores = ncores, 
                     parallel_workers = parallel_workers,
                     tune_length = tune_length, 
                     tune_stop = tune_stop, 
                     tune_time = tune_time, 
                     metric = metric, 
                     label = label, 
                     output = output, 
                     feature_type = feature_type, 
                     shap = shap)
      
    }
  }
}

pass_to_dietML <- function(train, test, metadata, model, program, seed, random_effects, folds, cv_repeats, ncores, parallel_workers, tune_length, tune_stop, tune_time, metric, label, output, feature_type, shap) {
  
  ## check for outdir and make if not there
  if (dir.exists(paste0(output, "/ml_analysis")) != TRUE) {
    dir.create(path = paste0(output, "/ml_analysis"))
  }
  
  
  ## check for label
  if ("feature_of_interest" %in% colnames(train) == FALSE & "feature_of_interest" %in% colnames(test) == FALSE) {
    stop(paste0("label not found in training AND testing data"))
  } 
  
  ## check if classification was mis-specified
  if (feature_type == "factor") {
    type <- "classification"
    if(length(levels(as.factor(metadata$feature_of_interest))) > 9)
      stop("You are trying to predict 10 or more classes. That is a bit much. Did you mean to do regression?")
  } else {
    type <- "regression"
  }
  
  ## run null model first, results_df is needed for the actually runs
  results_df <- run_null_model(train = train, test = test, seed = seed, type = type, random_effects = random_effects, output = output)
  
  ## if specified, run random forest
  if (model == "rf") {
    shap_inputs <- run_dietML_ranger(
      train = train,
      test = test,
      seed = seed,
      random_effects = random_effects,
      folds = folds,
      cv_repeats = cv_repeats,
      ncores = ncores,
      parallel_workers = parallel_workers,
      tune_length = tune_length,
      tune_stop = tune_stop,
      tune_time = tune_time,
      metric = metric,
      label = label,
      model = model,
      program = program,
      output = output,
      type = type,
      null_results = results_df
    )
  }
  if (model == "enet") {
    shap_inputs <- run_dietML_enet(
      train = train,
      test = test,
      seed = seed,
      random_effects = random_effects,
      folds = folds,
      cv_repeats = cv_repeats,
      ncores = ncores,
      parallel_workers = parallel_workers,
      tune_length = tune_length,
      tune_stop = tune_stop,
      tune_time = tune_time,
      metric = metric,
      label = label,
      model = model,
      program = program,
      output = output,
      type = type,
      null_results = results_df
    )
  }
  if (shap) {
    shap_analysis(label = label, 
                  output = output, 
                  model = model, 
                  filename = paste0(program, "_", seed), 
                  shap_inputs = shap_inputs,
                  train = train,
                  test = test,
                  type = type,
                  parallel_workers = parallel_workers)
  }
}

run_dietML_ranger <- function(train, test, seed, random_effects, folds, cv_repeats, parallel_workers, ncores, tune_length, tune_stop, tune_time, metric, label, model, program, output, type, null_results) {
  
  ## check and make sure results_df has been created from run_null_model,
  ## needed for this function downstream
  if (!exists("null_results")) {
    stop("Null model was not run. Cannot complete model evaluation.")
  }
  
  ## set total cores
  total_cores <- (as.numeric(ncores) * as.numeric(parallel_workers))
  
  ## remove individual and train if random effects
  if (random_effects) {
    train <- train %>% dplyr::select(., -dplyr::any_of(c("individual", "time")))
    test <- train %>% dplyr::select(., -dplyr::any_of(c("individual", "time")))
  }
  
  split_from_data_frame <- make_splits(
    x = train,
    assessment = test
  )
  
  ## set resampling scheme
  folds <- rsample::vfold_cv(train, v = as.numeric(folds), strata = feature_of_interest, repeats = cv_repeats)
  
  ## recipe
  
  ## specify recipe (this is like the pre-process work)
  diet_ml_recipe <- recipes::recipe(feature_of_interest ~ ., data = train) %>% 
    recipes::update_role(tidyr::any_of("subject_id"), new_role = "ID") %>% 
    recipes::step_zv(all_predictors()) %>%
    recipes::step_novel(all_nominal_predictors()) %>%
    recipes::step_dummy(recipes::all_nominal_predictors())
  
  ## ML engine
  
  ## specify ML model and engine 
  if (as.numeric(tune_time) == 0) {
    initial_mod <- parsnip::rand_forest(mode = type) %>%
      parsnip::set_engine("ranger", 
                          num.threads = total_cores,
                          importance = "none")
    
    initial_mod %>% parsnip::translate()
    
  } else {
    ## specify ML model and engine 
    initial_mod <- parsnip::rand_forest(mode = type, 
                                        mtry = tune(),
                                        trees = tune(),
                                        min_n = tune()) %>%
      parsnip::set_engine("ranger", 
                          num.threads = as.numeric(ncores),
                          importance = "none")
    
    initial_mod %>% parsnip::translate()
  }
  
  ## workflow ====================================================================
  
  ## define workflow
  dietML_wflow <- 
    workflows::workflow() %>% 
    workflows::add_model(initial_mod) %>% 
    workflows::add_recipe(diet_ml_recipe, blueprint = hardhat::default_recipe_blueprint(allow_novel_levels = TRUE))  
  #print(dietML_wflow)
  
  ## hyperparameters =============================================================
  
  if (as.numeric(tune_time) == 0) {
    no_tune_model <- parsnip::fit(dietML_wflow, train)
    ## create the last model based on best parameters
    last_best_mod <- 
      parsnip::rand_forest(mtry = no_tune_model$fit[[2]]$fit$mtry, min_n = no_tune_model$fit[[2]]$fit$min.node.size, trees = no_tune_model$fit[[2]]$fit$num.trees) %>% 
      parsnip::set_engine("ranger", num.threads = as.numeric(total_cores), importance = "none") %>% 
      parsnip::set_mode(type)
    
    ## update workflow with best model
    best_tidy_workflow <- 
      dietML_wflow %>% 
      workflows::update_model(last_best_mod)
    
  } else {
    ## define the hyper parameter set
    dietML_param_set <- parsnip::extract_parameter_set_dials(dietML_wflow)
    
    dietML_param_set <- 
      dietML_param_set %>% 
      # Pick an upper bound for mtry: 
      recipes::update(mtry = mtry(range(1, ncol(train %>% dplyr::select(., dplyr::any_of(c("feature_of_interest", "subject_id")))))))
    
    
    ## set up parallel jobs ========================================================
    ## remove any doParallel job setups that may have
    ## unneccessarily hung around
    unregister_dopar()
    
    ## register parallel cluster
    cl <- parallel::makePSOCKcluster(as.numeric(parallel_workers))
    doParallel::registerDoParallel(cl)
    
    ## set up hyper parameter search
    if (type == "classification") {
      
      search_res <-
        dietML_wflow %>% 
        tune::tune_bayes(
          resamples = folds,
          # To use non-default parameter ranges
          param_info = dietML_param_set,
          # Generate five at semi-random to start
          initial = 5,
          iter = tune_length,
          # How to measure performance?
          metrics = yardstick::metric_set(bal_accuracy, roc_auc, accuracy, kap, f_meas),
          control = tune::control_bayes(no_improve = as.numeric(tune_stop),
                                        uncertain = 5,
                                        verbose = FALSE,
                                        parallel_over = "resamples",
                                        time_limit = as.numeric(tune_time),
                                        seed = as.numeric(seed),
                                        save_pred = FALSE)
        )
      
    } else if (type == "regression") {
      
      search_res <-
        dietML_wflow %>% 
        tune::tune_bayes(
          resamples = folds,
          # To use non-default parameter ranges
          param_info = dietML_param_set,
          # Generate five at semi-random to start
          initial = 5,
          iter = tune_length,
          # How to measure performance?
          metrics = yardstick::metric_set(mae, rmse, rsq, ccc),
          control = tune::control_bayes(no_improve = as.numeric(tune_stop),
                                        uncertain = 5,
                                        verbose = FALSE,
                                        parallel_over = "resamples",
                                        time_limit = as.numeric(tune_time),
                                        seed = as.numeric(seed),
                                        save_pred = FALSE)
        )
    }
    
    ## stop parallel jobs
    parallel::stopCluster(cl)
    ## remove any doParallel job setups that may have
    ## unneccessarily hung around
    unregister_dopar()
    
    ## fit best model ============================================================
    
    ## get the best parameters from tuning
    best_mod <- 
      search_res %>% 
      tune::select_best(metric = metric)
    
    ## create the last model based on best parameters
    last_best_mod <- 
      parsnip::rand_forest(mtry = best_mod$mtry, min_n = best_mod$min_n, trees = best_mod$trees) %>% 
      parsnip::set_engine("ranger", num.threads = as.numeric(total_cores), importance = "none") %>% 
      parsnip::set_mode(type)
    
    ## update workflow with best model
    best_tidy_workflow <- 
      dietML_wflow %>% 
      workflows::update_model(last_best_mod)
    
    ## graphs ====================================================================
    
    hyperpar_perf_plot <- autoplot(search_res, type = "performance")
    ggplot2::ggsave(plot = hyperpar_perf_plot, filename = paste0(output, "/ml_analysis/", "training_performance.pdf"), width = 7, height = 2.5, units = "in")
    
    hyperpar_tested_plot <- autoplot(search_res, type = "parameters") + 
      labs(x = "Iterations", y = NULL)
    ggplot2::ggsave(plot = hyperpar_tested_plot, filename = paste0(output, "/ml_analysis/", "hyperpars_tested.pdf"), width = 7, height = 2.5, units = "in")
    
  }
  
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
  
  ## merge null results with trained results and write table
  null_results <- null_results %>% 
    dplyr::select(., -seed) %>% 
    summarise_all(., mean) %>% 
    t() %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = ".metric") %>% 
    dplyr::rename(., "null_model_avg" = 2)
  full_results <- merge(workflowsets::collect_metrics(final_res), null_results, by = ".metric", all = T)
  full_results$seed <- seed
  
  ## keep track of what program is being run for compete all levels
  full_results$program <- program
  
  
  ## write final results to file or append if file exists
  readr::write_csv(x = full_results, file = paste0(output, "/ml_analysis/ml_results.csv"), 
                   append = T, col_names = !file.exists(paste0(output, "/ml_analysis/ml_results.csv")))
  
  ## load up list for shap analysis
  shap_inputs <- list("split_from_data_frame" = split_from_data_frame, "diet_ml_recipe" = diet_ml_recipe, "best_tidy_workflow" = best_tidy_workflow)
  return(shap_inputs)
  
}

run_dietML_enet <- function(train, test, seed, random_effects, folds, cv_repeats, parallel_workers, ncores, tune_length, tune_stop, tune_time, metric, label, model, program, output, type, null_results) {
  
  ## check and make sure results_df has been created from run_null_model,
  ## needed for this function downstream
  if (!exists("null_results")) {
    stop("Null model was not run. Cannot complete model evaluation.")
  }
  
  ## set total cores
  total_cores <- (as.numeric(ncores) * as.numeric(parallel_workers))
  
  ## remove individual and train if random effects
  if (random_effects) {
    train <- train %>% dplyr::select(., -dplyr::any_of(c("individual", "time")))
    test <- train %>% dplyr::select(., -dplyr::any_of(c("individual", "time")))
  }
  
  split_from_data_frame <- make_splits(
    x = train,
    assessment = test
  )
  
  ## set resampling scheme
  folds <- rsample::vfold_cv(train, v = as.numeric(folds), strata = feature_of_interest, repeats = cv_repeats)
  
  ## recipe
  
  ## specify recipe (this is like the pre-process work)
  diet_ml_recipe <- recipes::recipe(feature_of_interest ~ ., data = train) %>% 
    recipes::update_role(tidyr::any_of("subject_id"), new_role = "ID") %>% 
    recipes::step_zv(all_predictors()) %>%
    recipes::step_novel(all_nominal_predictors()) %>%
    recipes::step_dummy(recipes::all_nominal_predictors())
  
  ## ML engine
  
  ## specify ML model and engine 
  if (as.numeric(tune_time) == 0) {
    # Define model with fixed penalty and mixture
    if (type == "classification") {
      initial_mod <- parsnip::logistic_reg(
        mode = "classification",
        penalty = double(1),
        mixture = 0.5
      ) %>%
        parsnip::set_engine("glmnet")
    } else {
      initial_mod <- parsnip::linear_reg(
        mode = "regression",
        penalty = double(1),
        mixture = 0.5
      ) %>%
        parsnip::set_engine("glmnet")
    }
  } else {
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
  } 
  
  initial_mod %>% parsnip::translate()
  
  ## define workflow
  dietML_wflow <- 
    workflows::workflow() %>% 
    workflows::add_model(initial_mod) %>% 
    workflows::add_recipe(diet_ml_recipe, blueprint = hardhat::default_recipe_blueprint(allow_novel_levels = TRUE))  
  #print(dietML_wflow)
  
  ## hyperparameters =============================================================
  if (as.numeric(tune_time) == 0) {
    best_tidy_workflow <- dietML_wflow 
  } 
  else {
    ## define the hyper parameter set
    dietML_param_set <- parsnip::extract_parameter_set_dials(dietML_wflow)
    
    ## set up parallel jobs ========================================================
    ## remove any doParallel job setups that may have
    ## unneccessarily hung around
    unregister_dopar()
    
    ## register parallel cluster
    cl <- parallel::makePSOCKcluster(as.numeric(parallel_workers))
    doParallel::registerDoParallel(cl)
    
    ## set up hyper parameter search
    if (type == "classification") {
      
      search_res <-
        dietML_wflow %>% 
        tune::tune_bayes(
          resamples = folds,
          # To use non-default parameter ranges
          param_info = dietML_param_set,
          # Generate five at semi-random to start
          initial = 5,
          iter = tune_length,
          # How to measure performance?
          metrics = yardstick::metric_set(bal_accuracy, roc_auc, accuracy, kap, f_meas),
          control = tune::control_bayes(no_improve = as.numeric(tune_stop),
                                        uncertain = 5,
                                        verbose = FALSE,
                                        parallel_over = "resamples",
                                        time_limit = as.numeric(tune_time),
                                        seed = as.numeric(seed),
                                        save_pred = FALSE)
        )
      
    } else if (type == "regression") {
      
      search_res <-
        dietML_wflow %>% 
        tune::tune_bayes(
          resamples = folds,
          # To use non-default parameter ranges
          param_info = dietML_param_set,
          # Generate five at semi-random to start
          initial = 5,
          iter = tune_length,
          # How to measure performance?
          metrics = yardstick::metric_set(mae, rmse, rsq, ccc),
          control = tune::control_bayes(no_improve = as.numeric(tune_stop),
                                        uncertain = 5,
                                        verbose = FALSE,
                                        parallel_over = "resamples",
                                        time_limit = as.numeric(tune_time),
                                        seed = as.numeric(seed),
                                        save_pred = FALSE)
        )
    }
    
    ## stop parallel jobs
    parallel::stopCluster(cl)
    ## remove any doParallel job setups that may have
    ## unneccessarily hung around
    unregister_dopar()
    
    ## fit best model ============================================================
    
    ## get the best parameters from tuning
    best_mod <- 
      search_res %>% 
      tune::select_best(metric = metric)
    
    ## create the last model based on best parameters
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
      dietML_wflow %>% 
      workflows::update_model(last_best_mod)
    
    ## graphs ====================================================================
    
    hyperpar_perf_plot <- autoplot(search_res, type = "performance")
    ggplot2::ggsave(plot = hyperpar_perf_plot, filename = paste0(output, "/ml_analysis/", "training_performance.pdf"), width = 7, height = 2.5, units = "in")
    
    hyperpar_tested_plot <- autoplot(search_res, type = "parameters") + 
      labs(x = "Iterations", y = NULL)
    ggplot2::ggsave(plot = hyperpar_tested_plot, filename = paste0(output, "/ml_analysis/", "hyperpars_tested.pdf"), width = 7, height = 2.5, units = "in")
    
  }
  
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
  
  ## merge null results with trained results and write table
  null_results <- null_results %>% 
    dplyr::select(., -seed) %>% 
    summarise_all(., mean) %>% 
    t() %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = ".metric") %>% 
    dplyr::rename(., "null_model_avg" = 2)
  full_results <- merge(workflowsets::collect_metrics(final_res), null_results, by = ".metric", all = T)
  full_results$seed <- seed
  
  ## keep track of what program is being run for compete all levels
  full_results$program <- program
  ## keep track of what model is being run for compete all levels
  full_results$model <- model
  
  ## write final results to file or append if file exists
  readr::write_csv(x = full_results, file = paste0(output, "/ml_analysis/ml_results.csv"), 
                   append = T, col_names = !file.exists(paste0(output, "/ml_analysis/ml_results.csv")))
  
  ## load up list for shap analysis
  shap_inputs <- list("split_from_data_frame" = split_from_data_frame, "diet_ml_recipe" = diet_ml_recipe, "best_tidy_workflow" = best_tidy_workflow)
  return(shap_inputs)
  
}

run_null_model <- function(train, test, seed, type, random_effects, output) {
  
  ## create results df
  if (type == "classification") {
    results_df <- data.frame(seed = numeric(), bal_accuracy = numeric(), f_meas = numeric(), accuracy = numeric(), stringsAsFactors = F)
  } else if (type == "regression") {
    results_df <- data.frame(seed = numeric(), mae = numeric(), rmse = numeric(), ccc = numeric(), stringsAsFactors = F)
  }
  
  ## interate over null model
  if (type == "classification") {
    df_loop_results <- data.frame(truth = character(), estimate = character(), stringsAsFactors = F)
  } else if (type == "regression") {
    df_loop_results <- data.frame(truth = numeric(), estimate = numeric(), stringsAsFactors = F)
  }
  
  ## remove individual and train if random effects
  if (random_effects) {
    train <- train %>% dplyr::select(., -dplyr::any_of(c("individual", "time")))
    test <- train %>% dplyr::select(., -dplyr::any_of(c("individual", "time")))
  }
  
  ## recipe
  
  ## specify recipe (this is like the pre-process work)
  diet_ml_recipe <- recipes::recipe(feature_of_interest ~ ., data = train) %>% 
    recipes::update_role(tidyr::any_of("subject_id"), new_role = "ID") %>% 
    recipes::step_zv(all_predictors()) %>%
    recipes::step_novel(all_nominal_predictors()) %>%
    recipes::step_dummy(recipes::all_nominal_predictors())
  
  
  ## ML engine
  
  ## specify ML model and engine 
  initial_mod <- null_model() %>% 
    set_engine("parsnip") %>% 
    set_mode(type) %>% 
    translate()
  
  ## workflow
  
  ## define workflow
  diet_ml_wflow <- 
    workflows::workflow() %>% 
    workflows::add_model(initial_mod) %>% 
    workflows::add_recipe(diet_ml_recipe, blueprint = hardhat::default_recipe_blueprint(allow_novel_levels = TRUE))  
  
  
  ## fit model
  
  ## fit to test data
  final_res <- parsnip::fit(diet_ml_wflow, test)
  
  df_loop_results <- add_row(df_loop_results, truth = test$feature_of_interest)
  df_loop_results$estimate <- final_res$fit$fit$fit$value
  
  if (type== "classification") {
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
                      seed = seed)
  } else if (type == "regression") {
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
                      seed = seed)
    
  }
  
  ## write table of results to file
  readr::write_csv(x = results_df, file =paste0(output, "/ml_analysis/dummy_model_results.csv"), 
                   append = T, col_names = !file.exists(paste0(output, "/ml_analysis/dummy_model_results.csv")))
  
  ## return results_df because that is what the other models need (ranger, enet)
  return(results_df)
  
}