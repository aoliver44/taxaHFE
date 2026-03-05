stacks_env <- new.env()
for (mod_type in c("enet", "ridge", "lasso", "svm", "rf", "xgboost", "mars")) {
  shap_inputs <- run_dietML(train = as.data.frame(train_data), 
                            test = as.data.frame(test_data), 
                            model = mod_type,
                            program = program, 
                            seed = opts$seed, 
                            random_effects = FALSE, 
                            nfolds = opts$folds, 
                            cv_repeats = opts$cv_repeats, 
                            ncores = opts$ncores, 
                            parallel_workers = opts$parallel_workers, 
                            tune_length = opts$tune_length, 
                            tune_stop = opts$tune_stop, 
                            tune_time = opts$tune_time, 
                            metric = opts$metric, 
                            label = opts$label, 
                            output = opts$output, 
                            feature_type = opts$feature_type, 
                            shap = opts$shap, 
                            cor_level = opts$cor_level, 
                            vif_threshold = opts$vif_threshold,
                            info_gain_n = opts$info_gain_n,
                            pct_loss = opts$pct_loss
  )
  #rm(shap_inputs)
}




stacks_attempt <- 
  stacks() |>
  add_candidates(stacks_env$enet_tune_res) |>
  add_candidates(stacks_env$mars_tune_res) |>
  add_candidates(stacks_env$rf_tune_res) |>
  add_candidates(stacks_env$svm_tune_res) |>
  add_candidates(stacks_env$xgboost_tune_res)
stacks_attempt_model_st <-
  stacks_attempt |>
  blend_predictions()
autoplot(stacks_attempt_model_st, type = "weights")
stacks_attempt_model_st_fit <-
  stacks_attempt_model_st |>
  fit_members()
test_data <- rsample::testing(tr_te_split)
test_data <- 
  bind_cols(test_data, predict(stacks_attempt_model_st_fit, test_data))
ggplot(test_data) +
  aes(
    x = feature_of_interest,
    y = .pred
  ) +
  geom_point() +
  coord_obs_pred()

yardstick::rsq_vec(test_data$feature_of_interest, test_data$.pred)
yardstick::mae_vec(test_data$feature_of_interest, test_data$.pred)

member_preds <- 
  test_data |>
  select(feature_of_interest) |>
  bind_cols(predict(stacks_attempt_model_st_fit, test_data, members = TRUE))
map(member_preds, rsq_vec, truth = member_preds$feature_of_interest) |>
  as_tibble()
