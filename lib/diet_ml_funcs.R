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
    no_tune_model <- parsnip::fit(diet_ml_workflow, split_from_data_frame$data[split_from_data_frame$in_id,])
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
  
  ## workflow ==================================================================
  
  ## define workflow
  diet_ml_workflow <- dietml_workflow(model_obj = initial_mod, recipe = diet_ml_recipe)
  
  ## hyperparameters =============================================================
  if (as.numeric(tune_time) == 0) {
    best_tidy_workflow <- diet_ml_workflow 
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
  ## specify regularizaton path for pure ridge regression because of theis issue
  ## https://github.com/tidymodels/parsnip/issues/431#issuecomment-782883848
  coef_path_values <- c(0, 10^seq(-6, 2, length.out = 100))
  
  ## specify ML model and engine 
  if (as.numeric(tune_time) == 0) {
    # Define model with fixed penalty and mixture
    if (type == "classification") {
      initial_mod <- parsnip::logistic_reg(
        mode = "classification",
        penalty = double(1),
        mixture = ifelse(model == "lasso", 1, 0)
      ) %>%
        {if (model == "lasso") parsnip::set_engine(., "glmnet") else parsnip::set_engine(., "glmnet", path_values = coef_path_values)}
    } else {
      initial_mod <- parsnip::linear_reg(
        mode = "regression",
        penalty = double(1),
        mixture = ifelse(model == "lasso", 1, 0)
      ) %>%
        {if (model == "lasso") parsnip::set_engine(., "glmnet") else parsnip::set_engine(., "glmnet", path_values = coef_path_values)}
    }
  } else {
    if (type == "classification") {
      initial_mod <- parsnip::logistic_reg(mode = "classification", 
                                           penalty = tune(),
                                           mixture = ifelse(model == "lasso", 1, 0)) %>%
        {if (model == "lasso") parsnip::set_engine(., "glmnet") else parsnip::set_engine(., "glmnet", path_values = coef_path_values)}
        
    } else {
      initial_mod <- parsnip::linear_reg(mode = "regression", 
                                         penalty = tune(),
                                         mixture = ifelse(model == "lasso", 1, 0)) %>%
        {if (model == "lasso") parsnip::set_engine(., "glmnet") else parsnip::set_engine(., "glmnet", path_values = coef_path_values)}
    }
  } 
  
  ## workflow ==================================================================
  
  ## define workflow
  diet_ml_workflow <- dietml_workflow(model_obj = initial_mod, recipe = diet_ml_recipe)
  
  ## hyperparameters =============================================================
  if (as.numeric(tune_time) == 0) {
    best_tidy_workflow <- diet_ml_workflow 
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
  diet_ml_workflow <- dietml_workflow(mode = initial_mod, recipe = diet_ml_recipe)
  ## fit model
  
  ## fit to test data
  final_res <- parsnip::fit(diet_ml_workflow, split_from_data_frame$data[split_from_data_frame$out_id,])
  
  null_estimates <- data.frame(estimate = final_res$fit$fit$fit$value, truth = split_from_data_frame$data[split_from_data_frame$out_id,]$feature_of_interest)
  
  if (type== "classification") {
    ## for yardstick, the estimate must have the same number of levels as
    ## the truth, even though the estimate will take on only one value
    null_estimates$estimate <- factor(x = null_estimates$estimate, levels = levels(as.factor(null_estimates$truth)))
    results_df <- results_df %>% 
      tibble::add_row(., bal_accuracy = 
                        yardstick::bal_accuracy_vec(truth = as.factor(null_estimates$truth), 
                                                    estimate = as.factor(null_estimates$estimate), 
                                                    data = null_estimates), 
                      accuracy = 
                        yardstick::accuracy_vec(truth = as.factor(null_estimates$truth), 
                                                estimate = as.factor(null_estimates$estimate), 
                                                data = null_estimates), 
                      f_meas = 
                        yardstick::f_meas_vec(truth = as.factor(null_estimates$truth), 
                                              estimate = as.factor(null_estimates$estimate), 
                                              data = null_estimates),
                      seed = seed)
  } else if (type == "regression") {
    results_df <- results_df %>% 
      tibble::add_row(., mae = 
                        yardstick::mae_vec(truth = null_estimates$truth, 
                                           estimate = null_estimates$estimate, 
                                           data = null_estimates), 
                      rmse = 
                        yardstick::rmse_vec(truth = null_estimates$truth, 
                                            estimate = null_estimates$estimate, 
                                            data = null_estimates),
                      ccc = yardstick::ccc_vec(truth = null_estimates$truth, 
                                               estimate = null_estimates$estimate, 
                                               data = null_estimates),
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
    test <- train %>% dplyr::select(., -dplyr::any_of(c("individual", "time")))
  }
  
  split_from_data_frame <- rsample::make_splits(x = train, assessment = test)
  return(split_from_data_frame)
}

set_cv_strategy <- function(split_from_data_frame, folds, feature_of_interest, cv_repeats) {
  train <- split_from_data_frame$data[split_from_data_frame$in_id,]
  ## set resampling scheme
  cv_folds <- rsample::vfold_cv(train, v = as.numeric(folds), strata = feature_of_interest, repeats = cv_repeats)
  
  ## log CV strategy
  logger::log_info("Stratified (across the response) cross validation strategy set, using {as.numeric(folds)} and repeating {cv_repeats}x time(s).")
  return(cv_folds)
}

dietml_recipe <- function(split_from_data_frame, cor_level, info_gain_n, type, ncores, model) {
  
  train <- split_from_data_frame$data[split_from_data_frame$in_id,]
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
  
  train <- split_from_data_frame$data[split_from_data_frame$in_id,]
  ## define the hyper parameter set
  dietML_param_set <- parsnip::extract_parameter_set_dials(diet_ml_workflow)
  
  ## make sure mtry is not ever more than the predictors we have
  if (model == "rf") {
    dietML_param_set <- 
      dietML_param_set %>% 
      # Pick an upper bound for mtry: 
      recipes::update(mtry = mtry(range(1, ncol(train %>% dplyr::select(., dplyr::any_of(c("feature_of_interest", "subject_id")))))))
  }
  
  ## make sure the penalty is a wide enough space, else some metrics like MAE
  ## can be the same, resulting in an error. Also we'll ask the tuner to 
  ## intialize more models to help prevent this.
  if (model %in% c("ridge", "lasso", "enet")) {
    dietML_param_set <- 
      dietML_param_set %>% 
      # widen the penalty search space, help prevent MAE from locking into zero variance: 
      recipes::update(penalty = penalty(range = c(-8, 3), trans = scales::transform_log10())) %>%
      {if (model == "enet") recipes::update(., mixture = mixture(range(0.1, 0.9))) else .}
    
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
                                      time_limit = as.numeric(tune_time),
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
  
  ## create the last model based on best parameters: Random Forest
  if (model == "rf") {
    last_best_mod <- 
      parsnip::rand_forest(mtry = best_mod$mtry, min_n = best_mod$min_n, trees = best_mod$trees) %>% 
      parsnip::set_engine("ranger", num.threads = as.numeric(ncores), importance = "none") %>% 
      parsnip::set_mode(type)
  } 
  
  ## create the last model based on best parameters: Elastic Net
  if (model == "enet") {
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
  }
  
  ## create the last model based on best parameters: Ridge or Lasso
  ## specify regularizaton path for pure ridge regression because of theis issue
  ## https://github.com/tidymodels/parsnip/issues/431#issuecomment-782883848
  coef_path_values <- c(0, 10^seq(-6, 2, length.out = 100))
  
  if (model %in% c("ridge", "lasso")) {
    ## create the last model based on best parameters
    if (type == "classification" && model == "lasso") {
      last_best_mod <- parsnip::logistic_reg(mode = "classification", penalty = best_mod$penalty, mixture = 1) %>% 
        parsnip::set_engine("glmnet") %>% 
        parsnip::set_mode(type)
    } else if (type == "classification" && model == "ridge") {
      last_best_mod <- parsnip::logistic_reg(mode = "classification", penalty = best_mod$penalty, mixture = 0) %>% 
        parsnip::set_engine("glmnet") %>% 
        parsnip::set_mode(type)
    } else if (type == "regression" && model == "lasso") {
      last_best_mod <- parsnip::linear_reg(mode = "regression", penalty = best_mod$penalty, mixture = 1) %>% 
        parsnip::set_engine("glmnet") %>% 
        parsnip::set_mode(type)
    } else if (type == "regression" && model == "ridge") {
      last_best_mod <- parsnip::linear_reg(mode = "regression", penalty = best_mod$penalty, mixture = 0) %>% 
        parsnip::set_engine("glmnet") %>% 
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
  
  ## return shap inputs list and results. The recipe and workflow are all 
  ## inside the final_res object. This is important to use because the preprocessing 
  ## steps were estimated on the training data and the model trained on the training data
  shap_inputs <- list("split_from_data_frame" = split_from_data_frame, "final_res" = final_res, "full_results" = full_results)
  return(shap_inputs)
  
}
