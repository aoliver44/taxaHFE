## DIETML FUNCTIONS

## libraries  ==================================================================
source("lib/requirements.R")

## helper functions ============================================================

## Negate function ("not in"):
`%!in%` <- Negate(`%in%`)

## suppress warnings
options(warn = -1)

## run dietML  ========================================

run_dietML <- function(train, test, model, program, seed, 
                       random_effects, folds, cv_repeats, ncores, 
                       parallel_workers, tune_length, tune_stop, tune_time, 
                       metric, label, output, feature_type, shap, cor_level, 
                       info_gain_n) {
  
  ## check for outdir and make if not there
  if (!dir.exists(paste0(output, "/ml_analysis"))) {
    dir.create(path = paste0(output, "/ml_analysis"))
  }

  ## check for label
  if ("feature_of_interest" %in% colnames(train) == FALSE & "feature_of_interest" %in% colnames(test) == FALSE) {
    logger::log_fatal("label not found in training AND testing data")
    stop()
  } 
  
  ## combine train and test data into a data split object
  split_from_data_frame <- create_data_split_obj(train = train, test = test, 
                                                 random_effects = random_effects)
  
  ## check if classification was mis-specified
  if (feature_type == "factor") {
    type <- "classification"
    if (length(levels(as.factor(split_from_data_frame$data$feature_of_interest))) > 9) {
      logger::log_fatal("You are trying to predict 10 or more classes. That is a bit much. Did you mean to do regression?")
      stop() 
      } 
    } else {
    type <- "regression"
  }
  
  ## run null model first, results_df is needed for the actually runs
  null_results <- run_null_model(split_from_data_frame = split_from_data_frame, seed = seed, 
                                 type = type, output = output, cv_repeats = cv_repeats, 
                                 feature_of_interest = "feature_of_interest", folds = folds, cor_level = cor_level, 
                                 info_gain_n = info_gain_n, ncores = ncores, model = model)
  
  ## if specified, run random forest
  if (model == "rf") {
    shap_inputs <- run_dietML_ranger(split_from_data_frame = split_from_data_frame, 
                                     seed = seed, folds = folds, cv_repeats = cv_repeats, 
                                     parallel_workers = parallel_workers, ncores = ncores, 
                                     tune_length = tune_length, tune_stop = tune_stop, 
                                     tune_time = tune_time, metric = metric, 
                                     model = model, program = program, output = output, 
                                     feature_of_interest = "feature_of_interest",
                                     type = type, null_results = null_results,
                                     cor_level = cor_level, info_gain_n = info_gain_n
                                     )
  }
  if (model == "enet") {
    shap_inputs <- run_dietML_enet(split_from_data_frame = split_from_data_frame, 
                                   seed = seed, folds = folds, cv_repeats = cv_repeats, 
                                   parallel_workers = parallel_workers, ncores = ncores, 
                                   tune_length = tune_length, tune_stop = tune_stop, 
                                   tune_time = tune_time, metric = metric, 
                                   model = model, program = program, output = output, 
                                   feature_of_interest = "feature_of_interest",
                                   type = type, null_results = null_results,
                                   cor_level = cor_level, info_gain_n = info_gain_n
    )
  }
  if (model %in% c("ridge", "lasso")) {
    shap_inputs <- run_dietML_ridge_lasso(split_from_data_frame = split_from_data_frame, 
                                   seed = seed, folds = folds, cv_repeats = cv_repeats, 
                                   parallel_workers = parallel_workers, ncores = ncores, 
                                   tune_length = tune_length, tune_stop = tune_stop, 
                                   tune_time = tune_time, metric = metric, 
                                   model = model, program = program, output = output, 
                                   feature_of_interest = "feature_of_interest",
                                   type = type, null_results = null_results,
                                   cor_level = cor_level, info_gain_n = info_gain_n
    )
  }
  return(shap_inputs)
}

run_dietML_ranger <- function(split_from_data_frame, seed, folds, cv_repeats, 
                              parallel_workers, ncores, tune_length, tune_stop, 
                              tune_time, metric, feature_of_interest, model, program, 
                              output, type, null_results, cor_level, info_gain_n) {
  
  ## log start of RF function
  logger::log_info("{model} model started...")
  
  ## set resampling scheme
  folds <- set_cv_strategy(split_from_data_frame = split_from_data_frame, 
                           folds = folds, feature_of_interest = feature_of_interest, 
                           cv_repeats = cv_repeats)
    
  ## recipe
  diet_ml_recipe <- dietml_recipe(split_from_data_frame = split_from_data_frame, 
                                  cor_level = cor_level, info_gain_n = info_gain_n, 
                                  type = type, ncores = ncores, model = model)
  
  ## Random Forest ML engine
  ## specify ML model and engine 
  ## if no HP tuning, set initial RF model, which will take defaults
  if (as.numeric(tune_time) == 0) {
    initial_mod <- parsnip::rand_forest(mode = type) %>%
      parsnip::set_engine("ranger", 
                          num.threads = as.numeric(ncores),
                          importance = "none")
    ## else if HP tuning time, set parameters to tune
  } else {
    initial_mod <- parsnip::rand_forest(mode = type, 
                                        mtry = tune(),
                                        trees = tune(),
                                        min_n = tune()) %>%
      parsnip::set_engine("ranger", 
                          num.threads = as.numeric(ncores),
                          importance = "none")
  }
  
  ## workflow ====================================================================
  
  ## define workflow
  diet_ml_workflow <- dietml_workflow(model_obj = initial_mod, recipe = diet_ml_recipe)
  
  ## hyperparameters =============================================================
  
  if (as.numeric(tune_time) == 0) {
    no_tune_model <- parsnip::fit(diet_ml_workflow, rsample::training(split_from_data_frame))
    ## create the last model based on best parameters
    last_best_mod <- 
      parsnip::rand_forest(mtry = no_tune_model$fit[[2]]$fit$mtry, min_n = no_tune_model$fit[[2]]$fit$min.node.size, trees = no_tune_model$fit[[2]]$fit$num.trees) %>% 
      parsnip::set_engine("ranger", num.threads = as.numeric(ncores), importance = "none") %>% 
      parsnip::set_mode(type)
    
    ## update workflow with best model
    best_tidy_workflow <- 
      diet_ml_workflow %>% 
      workflows::update_model(last_best_mod)
    
  } else {
    best_tidy_workflow <- dietml_hp_tune(split_from_data_frame = split_from_data_frame, 
                                         diet_ml_workflow = diet_ml_workflow, model = model, 
                                         parallel_workers = parallel_workers, 
                                         folds = folds, type = type, tune_time = tune_time, 
                                         seed = seed, tune_stop = tune_stop, metric = metric, 
                                         ncores = ncores, output = output, tune_length = tune_length)
  }
  
  ## write dietml outputs
  shap_inputs <- write_dietml_outputs(type = type, best_tidy_workflow = best_tidy_workflow, 
                                       split_from_data_frame = split_from_data_frame,
                                       seed = seed, null_results = null_results, 
                                       program = program, output = output)
  
  ## log end of rf model
  logger::log_info("{model} model finished!")
  
  ## return outputs
  return(shap_inputs)
  
}

run_dietML_enet <- function(split_from_data_frame, seed, folds, cv_repeats, 
                            parallel_workers, ncores, tune_length, tune_stop, 
                            tune_time, metric, feature_of_interest, model, program, 
                            output, type, null_results, cor_level, info_gain_n) {
  
  ## log start of ENET function
  logger::log_info("{model} model started...")
  
  ## set resampling scheme
  folds <- set_cv_strategy(split_from_data_frame = split_from_data_frame, 
                           folds = folds, feature_of_interest = feature_of_interest, 
                           cv_repeats = cv_repeats)
  
  ## recipe
  diet_ml_recipe <- dietml_recipe(split_from_data_frame = split_from_data_frame, 
                                  cor_level = cor_level, info_gain_n = info_gain_n, 
                                  type = type, ncores = ncores, model = model)
  
  ## ML engine
  ## specify ML model and engine 
  ## if tune_time = 0, only lightly tune penalty, which is mandatory.
    # Define model with fixed penalty and mixture
    if (type == "classification" && length(levels(as.factor(split_from_data_frame$data$feature_of_interest))) == 2) {
      initial_mod <- parsnip::logistic_reg(
        mode = "classification",
        penalty = tune(),
        mixture = ifelse(tune_time == 0, 0.5, tune())
      ) %>%
        parsnip::set_engine("glmnet") 
    } else if (type == "classification" && length(levels(as.factor(split_from_data_frame$data$feature_of_interest))) > 2) {
      initial_mod <- parsnip::multinom_reg(
        mode = "classification",
        penalty = tune(),
        mixture = ifelse(tune_time == 0, 0.5, tune())
      ) %>%
        parsnip::set_engine("glmnet") 
    } else if (type == "regression") {
      initial_mod <- parsnip::linear_reg(
        mode = "regression",
        penalty = tune(),
        mixture = ifelse(tune_time == 0, 0.5, tune())
      ) %>%
        parsnip::set_engine("glmnet")
    }
  
  ## workflow ==================================================================
  
  ## define workflow
  diet_ml_workflow <- dietml_workflow(model_obj = initial_mod, recipe = diet_ml_recipe)
  
  ## hyperparameters =============================================================
  if (as.numeric(tune_time) == 0) {
    logger::log_warn("You specified tune_time = 0, but penalized regression requires selecting a penalty. 
                     We do a fast, light cross-validation to choose a reasonable value. However, if your 
                     dataset is large, this may take a little time.")
    best_tidy_workflow <- dietml_hp_tune(split_from_data_frame = split_from_data_frame, 
                                         diet_ml_workflow = diet_ml_workflow, model = model, 
                                         parallel_workers = parallel_workers, 
                                         folds = folds, type = type, tune_time = 0, 
                                         seed = seed, tune_stop = tune_stop, metric = metric, 
                                         ncores = ncores, output = output, tune_length = 10)
  } 
  else {
    best_tidy_workflow <- dietml_hp_tune(split_from_data_frame = split_from_data_frame, 
                                         diet_ml_workflow = diet_ml_workflow, model = model, 
                                         parallel_workers = parallel_workers, 
                                         folds = folds, type = type, tune_time = tune_time, 
                                         seed = seed, tune_stop = tune_stop, metric = metric, 
                                         ncores = ncores, output = output, tune_length = tune_length)
  }
    
  ## write dietml outputs
  shap_inputs <- write_dietml_outputs(type = type, best_tidy_workflow = best_tidy_workflow, 
                                       split_from_data_frame = split_from_data_frame, 
                                       seed = seed, null_results = null_results, 
                                       program = program, output = output)
  
  ## log end of enet model
  logger::log_info("{model} model finished!")
  
  ## return outputs
  return(shap_inputs)
  
}

run_dietML_ridge_lasso <- function(split_from_data_frame, seed, folds, cv_repeats, 
                                   parallel_workers, ncores, tune_length, tune_stop, 
                                   tune_time, metric, feature_of_interest, model, program, 
                                   output, type, null_results, cor_level, info_gain_n) {
  
  ## log start of ridge, lasso function
  logger::log_info("{model} model started...")
  
  ## set resampling scheme
  folds <- set_cv_strategy(split_from_data_frame = split_from_data_frame, 
                           folds = folds, feature_of_interest = feature_of_interest, 
                           cv_repeats = cv_repeats)
  
  ## recipe
  diet_ml_recipe <- dietml_recipe(split_from_data_frame = split_from_data_frame, 
                                  cor_level = cor_level, info_gain_n = info_gain_n, 
                                  type = type, ncores = ncores, model = model)
  
  ## ML engine
  ## specify ML model and engine 
  if (type == "classification" && length(levels(as.factor(split_from_data_frame$data$feature_of_interest))) == 2) {
      initial_mod <- parsnip::logistic_reg(mode = "classification", 
                                           penalty = tune(),
                                           mixture = ifelse(model == "lasso", 1, 0)) %>%
        parsnip::set_engine("glmnet", standardize = FALSE)
    } else if (type == "classification" && length(levels(as.factor(split_from_data_frame$data$feature_of_interest))) > 2) {
      initial_mod <- parsnip::multinom_reg(
        mode = "classification",
        penalty = tune(),
        mixture = ifelse(model == "lasso", 1, 0)
      ) %>%
        parsnip::set_engine("glmnet") 
    } else if (type == "regression") {
      initial_mod <- parsnip::linear_reg(mode = "regression", 
                                         penalty = tune(),
                                         mixture = ifelse(model == "lasso", 1, 0)) %>%
        parsnip::set_engine("glmnet", standardize = FALSE)
    }
  
  ## workflow ==================================================================
  
  ## define workflow
  diet_ml_workflow <- dietml_workflow(model_obj = initial_mod, recipe = diet_ml_recipe)
  
  ## hyperparameters =============================================================
  if (as.numeric(tune_time) == 0) {
    logger::log_warn("You specified tune_time = 0, but penalized regression requires selecting a penalty. 
                     We do a fast, light cross-validation to choose a reasonable value. However, if your 
                     dataset is large, this may take a little time.")
    best_tidy_workflow <- dietml_hp_tune(split_from_data_frame = split_from_data_frame, 
                                         diet_ml_workflow = diet_ml_workflow, model = model, 
                                         parallel_workers = parallel_workers, 
                                         folds = folds, type = type, tune_time = 0, 
                                         seed = seed, tune_stop = tune_stop, metric = metric, 
                                         ncores = ncores, output = output, tune_length = 10)
  } 
  else {
    best_tidy_workflow <- dietml_hp_tune(split_from_data_frame = split_from_data_frame, 
                                         diet_ml_workflow = diet_ml_workflow, model = model, 
                                         parallel_workers = parallel_workers, 
                                         folds = folds, type = type, tune_time = tune_time, 
                                         seed = seed, tune_stop = tune_stop, metric = metric, 
                                         ncores = ncores, output = output, tune_length = tune_length)
  }
  
  ## write dietml outputs
  shap_inputs <- write_dietml_outputs(type = type, best_tidy_workflow = best_tidy_workflow, 
                                       split_from_data_frame = split_from_data_frame, 
                                       seed = seed, null_results = null_results, 
                                       program = program, output = output)
  
  ## log end of ridge/lasso model
  logger::log_info("{model} model finished!")
  
  ## return outputs
  return(shap_inputs)
  
}
run_null_model <- function(split_from_data_frame, seed, type, output, cv_repeats, feature_of_interest, folds, cor_level, info_gain_n, ncores, model) {
  
  ## log start of null model
  logger::log_info("null model started...")
  
  ## create results df
  if (type == "classification") {
    results_df <- data.frame(seed = numeric(), bal_accuracy = numeric(), f_meas = numeric(), accuracy = numeric(), stringsAsFactors = F)
  } else if (type == "regression") {
    results_df <- data.frame(seed = numeric(), mae = numeric(), rmse = numeric(), ccc = numeric(), stringsAsFactors = F)
  }
  
  ## set resampling scheme
  folds <- set_cv_strategy(split_from_data_frame = split_from_data_frame, 
                           folds = folds, feature_of_interest = feature_of_interest, 
                           cv_repeats = cv_repeats)
  
  ## recipe
  diet_ml_recipe <- dietml_recipe(split_from_data_frame = split_from_data_frame, 
                                  cor_level = cor_level, info_gain_n = info_gain_n, 
                                  type = type, ncores = ncores, model = model)
  
  ## specify ML model and engine 
  initial_mod <- null_model() %>% 
    set_engine("parsnip") %>% 
    set_mode(type) %>% 
    translate()
  
  ## workflow ==================================================================
  
  ## define workflow
  diet_ml_workflow <- dietml_workflow(model_obj = initial_mod, recipe = diet_ml_recipe)
  
  ## fit to test data using last_fit, same as the other models
  if (type == "classification") {
    final_res <- tune::last_fit(diet_ml_workflow, split_from_data_frame,
                                metrics = yardstick::metric_set(bal_accuracy, 
                                                                roc_auc, accuracy, 
                                                                kap, f_meas))
  } else if (type == "regression") {
    final_res <- tune::last_fit(diet_ml_workflow, split_from_data_frame,
                                metrics = yardstick::metric_set(mae, rmse, rsq, 
                                                                ccc))
  }
  
  ## extract predictions on test set
  null_estimates <- workflowsets::collect_predictions(final_res)
  
  if (type == "classification") {
    null_estimates$feature_of_interest <- factor(null_estimates$feature_of_interest)
    null_estimates$.pred_class <- factor(null_estimates$.pred_class, levels = levels(null_estimates$feature_of_interest))
    results_df <- results_df %>% 
      tibble::add_row(., bal_accuracy = 
                        yardstick::bal_accuracy_vec(truth = null_estimates$feature_of_interest, 
                                                    estimate = null_estimates$.pred_class), 
                      accuracy = 
                        yardstick::accuracy_vec(truth = null_estimates$feature_of_interest, 
                                                estimate = null_estimates$.pred_class), 
                      f_meas = 
                        yardstick::f_meas_vec(truth = null_estimates$feature_of_interest, 
                                              estimate = null_estimates$.pred_class),
                      seed = seed)
  } else if (type == "regression") {
    results_df <- results_df %>% 
      tibble::add_row(., mae = 
                        yardstick::mae_vec(truth = null_estimates$feature_of_interest, 
                                           estimate = null_estimates$.pred), 
                      rmse = 
                        yardstick::rmse_vec(truth = null_estimates$feature_of_interest, 
                                            estimate = null_estimates$.pred),
                      ccc = yardstick::ccc_vec(truth = null_estimates$feature_of_interest, 
                                               estimate = null_estimates$.pred),
                      seed = seed)
  }
  
  ## write table of results to file
  readr::write_csv(x = results_df, file = paste0(output, "/ml_analysis/dummy_model_results.csv"), 
                   append = T, col_names = !file.exists(paste0(output, "/ml_analysis/dummy_model_results.csv")))
  
  ## log end of null model
  logger::log_info("null model finished!")
  
  ## return results_df because that is what the other models need (ranger, enet)
  return(results_df)
}

create_data_split_obj <- function(train, test, random_effects) {

  ## remove individual and train if random effects
  if (random_effects) {
    train <- train %>% dplyr::select(., -dplyr::any_of(c("individual", "time")))
    test <- test %>% dplyr::select(., -dplyr::any_of(c("individual", "time")))
  }
  
  split_from_data_frame <- rsample::make_splits(x = train, assessment = test)
  return(split_from_data_frame)
}

set_cv_strategy <- function(split_from_data_frame, folds, feature_of_interest, cv_repeats) {
  train <- rsample::training(split_from_data_frame)
  ## set resampling scheme
  cv_folds <- rsample::vfold_cv(train, v = as.numeric(folds), strata = feature_of_interest, repeats = cv_repeats)
  
  ## log CV strategy
  logger::log_info("Stratified (across the response) cross validation strategy set, using {as.numeric(folds)} folds and repeating {cv_repeats}x time(s).")
  return(cv_folds)
}

dietml_recipe <- function(split_from_data_frame, cor_level, info_gain_n, type, ncores, model) {
  
  train <- rsample::training(split_from_data_frame)
  ## specify recipe (this is like the pre-process work)
  dietML_recipe <- recipes::recipe(feature_of_interest ~ ., data = train) %>% 
    recipes::update_role("subject_id", new_role = "ID") %>% 
    recipes::step_novel(recipes::all_nominal_predictors()) %>%
    recipes::step_dummy(recipes::all_nominal_predictors()) %>% 
    recipes::step_zv(recipes::all_predictors()) %>%
    {if (model %in% c("ridge", "lasso", "enet")) recipes::step_center(., recipes::all_numeric_predictors()) %>% 
        recipes::step_scale(., recipes::all_numeric_predictors()) else .} %>%
    {if (cor_level < 1) recipes::step_corr(., recipes::all_numeric_predictors(), threshold = cor_level, use = "everything") else .} %>%
    {if (info_gain_n > 0) colino::step_select_infgain(., recipes::all_predictors(), 
                                                          top_p = info_gain_n,
                                                          outcome = "feature_of_interest",
                                                          threads = ncores, scores = "tmp_scores") else .}
  
  ## idea - log intermediate file of what these steps do to the data
  
  return(dietML_recipe)
  
}

dietml_workflow <- function(model_obj, recipe) {
  ## define workflow
  dietML_wflow <- 
    workflows::workflow() %>% 
    workflows::add_model(model_obj) %>% 
    workflows::add_recipe(recipe, blueprint = hardhat::default_recipe_blueprint(allow_novel_levels = TRUE)) 
  
}

dietml_hp_tune <- function(diet_ml_workflow, model, parallel_workers, folds, type, tune_time, seed, tune_stop, metric, ncores, output, split_from_data_frame, tune_length) {
  
  ## for tuning, initiate the search space with 5 random models unless
  ## told otherwise
  n_inital_models = 5
  
  train <- rsample::training(split_from_data_frame)
  ## define the hyper parameter set
  dietML_param_set <- parsnip::extract_parameter_set_dials(diet_ml_workflow)
  
  ## make sure mtry is not ever more than the predictors we have
  if (model == "rf") {
    dietML_param_set <- 
      dietML_param_set %>% 
      # Pick an upper bound for mtry: 
      recipes::update(mtry = mtry(range(1, ncol(train %>% dplyr::select(., -dplyr::any_of(c("feature_of_interest", "subject_id"))))))) 
  }
  
  ## make sure the penalty is a wide enough space, else some metrics like MAE
  ## can be the same, resulting in an error. Also we'll ask the tuner to 
  ## intialize more models to help prevent this.
  if (model %in% c("ridge", "lasso", "enet")) {
    dietML_param_set <- 
      dietML_param_set %>% 
      # widen the penalty search space, help prevent MAE from locking into zero variance: 
      recipes::update(penalty = penalty(range = c(-8, 3), trans = scales::transform_log10())) %>%
      {if (model == "enet" && tune_time > 0) recipes::update(., mixture = mixture(range(0.1, 0.9))) else .}
    
    n_inital_models = 20
  }
  
  ## set up parallel jobs
  ## remove any doParallel job setups that may have
  ## unneccessarily hung around
  unregister_dopar()
  
  ## register parallel cluster
  cl <- parallel::makePSOCKcluster(as.numeric(parallel_workers))
  doParallel::registerDoParallel(cl)
  
  ## set up hyper parameter search
  if (type == "classification") {
    
    search_res <-
      diet_ml_workflow %>% 
      tune::tune_bayes(
        resamples = folds,
        # To use non-default parameter ranges
        param_info = dietML_param_set,
        # Generate five at semi-random to start
        initial = n_inital_models,
        iter = tune_length,
        # How to measure performance?
        metrics = yardstick::metric_set(bal_accuracy, roc_auc, accuracy, kap, f_meas),
        control = tune::control_bayes(no_improve = as.numeric(tune_stop),
                                      uncertain = 5,
                                      verbose = FALSE,
                                      parallel_over = "resamples",
                                      time_limit = ifelse(tune_time == 0 && model %in% c("enet", "ridge", "lasso"), 10e5, as.numeric(tune_time)), # fake super long time to allow 10 iterations to find penalty for tune_time = 0 (penalized regression)
                                      seed = as.numeric(seed),
                                      save_pred = FALSE)
      )
    
  } else if (type == "regression") {
    
    search_res <-
      diet_ml_workflow %>% 
      tune::tune_bayes(
        resamples = folds,
        # To use non-default parameter ranges
        param_info = dietML_param_set,
        # Generate five at semi-random to start
        initial = n_inital_models,
        iter = tune_length,
        # How to measure performance?
        metrics = yardstick::metric_set(mae, rmse, rsq, ccc),
        control = tune::control_bayes(no_improve = as.numeric(tune_stop),
                                      uncertain = 5,
                                      verbose = FALSE,
                                      parallel_over = "resamples",
                                      time_limit = ifelse(tune_time == 0 && model %in% c("enet", "ridge", "lasso"), 10e5, as.numeric(tune_time)), # fake super long time to allow 10 iterations to find penalty for tune_time = 0 (penalized regression)
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
  
  ## can also get a more regularized tree?
  # best_mod_other <- 
  #   search_res %>% 
  #   tune::select_by_one_std_err(metric = metric, desc(min_n), desc(mtry), desc(trees))
  
  ## log hyperparameters selected for best
  logger::log_info("Hyperparameters selected: ")
  for (n in seq(1:ncol(best_mod))) {
    logger::log_info(paste0(colnames(best_mod[n]), ": ", best_mod[1,n]))
  }
  
  ## create the last model based on best parameters: Random Forest
  if (model == "rf") {
    last_best_mod <- 
      parsnip::rand_forest(mtry = best_mod$mtry, min_n = best_mod$min_n, trees = best_mod$trees) %>% 
      parsnip::set_engine("ranger", num.threads = as.numeric(ncores), importance = "none") %>% 
      parsnip::set_mode(type)
  } 
  
  ## create the last model based on best parameters: Elastic Net
  ## if no tune time, make sure mixture is is set to 0.5
  if (model == "enet") {
    if (type == "classification" && length(levels(as.factor(split_from_data_frame$data$feature_of_interest))) == 2) {
      last_best_mod <- 
        parsnip::logistic_reg(mode = "classification", penalty = best_mod$penalty, mixture = ifelse(tune_time == 0, 0.5, best_mod$mixture)) %>% 
        parsnip::set_engine("glmnet", standardize = FALSE) %>% 
        parsnip::set_mode(type)
    } else if (type == "classification" && length(levels(as.factor(split_from_data_frame$data$feature_of_interest))) > 2) { 
      last_best_mod <- 
        parsnip::multinom_reg(mode = "classification", penalty = best_mod$penalty, mixture = ifelse(tune_time == 0, 0.5, best_mod$mixture)) %>% 
        parsnip::set_engine("glmnet", standardize = FALSE) %>% 
        parsnip::set_mode(type)
      } else if(type == "regression") {
      last_best_mod <- 
        parsnip::linear_reg(mode = "regression", penalty = best_mod$penalty, mixture = ifelse(tune_time == 0, 0.5, best_mod$mixture)) %>% 
        parsnip::set_engine("glmnet", standardize = FALSE) %>% 
        parsnip::set_mode(type)
    }
  }
  
  ## create the last model based on best parameters: Ridge or Lasso
  if (model %in% c("ridge", "lasso")) {
    ## create the last model based on best parameters
    if (type == "classification" && model == "lasso") {
      if (length(levels(as.factor(split_from_data_frame$data$feature_of_interest))) == 2) {
        last_best_mod <- parsnip::logistic_reg(mode = "classification", penalty = best_mod$penalty, mixture = 1) %>% 
          parsnip::set_engine("glmnet", standardize = FALSE) %>% 
          parsnip::set_mode(type) 
      } else if (length(levels(as.factor(split_from_data_frame$data$feature_of_interest))) > 2) {
        last_best_mod <- parsnip::multinom_reg(mode = "classification", penalty = best_mod$penalty, mixture = 1) %>% 
          parsnip::set_engine("glmnet", standardize = FALSE) %>% 
          parsnip::set_mode(type) 
      }
    } else if (type == "classification" && model == "ridge") {
      if (length(levels(as.factor(split_from_data_frame$data$feature_of_interest))) == 2) {
        last_best_mod <- parsnip::logistic_reg(mode = "classification", penalty = best_mod$penalty, mixture = 0) %>% 
          parsnip::set_engine("glmnet", standardize = FALSE) %>% 
          parsnip::set_mode(type) 
      } else if (length(levels(as.factor(split_from_data_frame$data$feature_of_interest))) > 2) {
        last_best_mod <- parsnip::multinom_reg(mode = "classification", penalty = best_mod$penalty, mixture = 0) %>% 
          parsnip::set_engine("glmnet", standardize = FALSE) %>% 
          parsnip::set_mode(type) 
      }
    } else if (type == "regression" && model == "lasso") {
      last_best_mod <- parsnip::linear_reg(mode = "regression", penalty = best_mod$penalty, mixture = 1) %>% 
        parsnip::set_engine("glmnet", standardize = FALSE) %>% 
        parsnip::set_mode(type)
    } else if (type == "regression" && model == "ridge") {
      last_best_mod <- parsnip::linear_reg(mode = "regression", penalty = best_mod$penalty, mixture = 0) %>% 
        parsnip::set_engine("glmnet", standardize = FALSE) %>% 
        parsnip::set_mode(type)
    }
  }
  
  ## update workflow with best model
  best_tidy_workflow <- 
    diet_ml_workflow %>% 
    workflows::update_model(last_best_mod)
  
  ## graphs ====================================================================
  
  hyperpar_perf_plot <- autoplot(search_res, type = "performance")
  ggplot2::ggsave(plot = hyperpar_perf_plot, filename = paste0(output, "/ml_analysis/", "training_performance.pdf"), width = 7, height = 2.5, units = "in")
  
  hyperpar_tested_plot <- autoplot(search_res, type = "parameters") + 
    labs(x = "Iterations", y = NULL)
  ggplot2::ggsave(plot = hyperpar_tested_plot, filename = paste0(output, "/ml_analysis/", "hyperpars_tested.pdf"), width = 7, height = 2.5, units = "in")
  
  return(best_tidy_workflow)
}

write_dietml_outputs <- function(type,  best_tidy_workflow, split_from_data_frame, seed, null_results, program, output) {
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
  
  ## calculate training scores and log them 
  all_predictions <- bind_rows(tune::augment(final_res$.workflow[[1]], new_data = rsample::training(split_from_data_frame)) %>% 
                                 dplyr::mutate(.model_input_type = "train"), 
                               tune::augment(final_res$.workflow[[1]], new_data = rsample::testing(split_from_data_frame)) %>% 
                                 dplyr::mutate(.model_input_type = "test"))
  if (type == "classification") {
    class_metrics <- yardstick::metric_set(bal_accuracy, accuracy, kap, f_meas)
    all_predictions <- all_predictions %>% dplyr::group_by(.model_input_type) %>% 
      dplyr::mutate(feature_of_interest = factor(feature_of_interest), .pred_class = factor(.pred_class)) %>% 
      class_metrics(truth = feature_of_interest, estimate = .pred_class)
  } else {
    all_predictions <- all_predictions %>% dplyr::group_by(.model_input_type) %>% yardstick::metrics(feature_of_interest, .pred)
  }
  
  logger::log_info("ML results (assessment on test data recorded in ml_results.csv):")
  for (n in seq(1:nrow(all_predictions))) {
    logger::log_info(paste0(all_predictions[n,2], " (", all_predictions[n,1], "): ", all_predictions[n,4]))
  }
  
  ## return shap inputs list and results. The recipe and workflow are all 
  ## inside the final_res object. This is important to use because the preprocessing 
  ## steps were estimated on the training data and the model trained on the training data
  shap_inputs <- list("split_from_data_frame" = split_from_data_frame, "final_res" = final_res, "full_results" = full_results)
  return(shap_inputs)
  
}
