## SHAP FUNCTIONS

## libraries  ==================================================================
source("lib/requirements.R")

## helper functions ============================================================

## Negate function ("not in"):
`%!in%` <- Negate(`%in%`)

## suppress warnings
options(warn = -1)

shap_analysis <- function(label, output, model, filename, shap_inputs, train, test, feature_type, parallel_workers) {
  
  # --- Setup ---
  shap_plot_env <- new.env()
  shap.error.occured <- FALSE
  error_message <- NULL
  output_dir <- paste0(output, "/ml_analysis")
  
  # --- Load SHAP inputs ---
  split_from_data_frame <- shap_inputs$split_from_data_frame
  final_workflow <- shap_inputs$final_res$.workflow[[1]]
  prepped_recipe <- workflows::extract_recipe(final_workflow)

  ## save some initial inputs to env, in case the below 
  ## shap analysis does not finish. Occasionaly it does not finish on
  ## the "test" dataset. Which is fine, i cant think of why that is used.
  ## But if it fails, we still want as much data returned as possible, so 
  ## that is why we return everything prior to returning the test shap data
  assign("split_from_data_frame", split_from_data_frame, envir = shap_plot_env)
  assign("final_workflow", final_workflow, envir = shap_plot_env)
  assign("label", label, envir = shap_plot_env)
  
  # --- Define prediction wrapper (pfun) ---
  pfun <- NULL
  pfun <- function(object, newdata) {
    
    # Binary classification
    if (feature_type == "factor" && length(levels(as.factor(split_from_data_frame$data$feature_of_interest))) == 2) {
      
      preds <- predict(object, new_data = newdata, type = "prob")
      pos_class <- levels(as.factor(split_from_data_frame$data$feature_of_interest))[1]
      logger::log_info(paste0("Positive class (pfun): ", pos_class))
      return(as.numeric(preds[[paste0(".pred_", pos_class)]]))
      
    } else if (feature_type == "numeric") {
      # Regression
      preds <- predict(object, new_data = newdata, type = "numeric")
      as.numeric(preds$.pred)
    }
  }
  
  if (is.null(pfun)) {
    message("Error: Could not define prediction function (pfun). Check model and type inputs.")
    shap.error.occured <- TRUE
  } else {
    # --- SHAP analysis block ---
    result <- tryCatch({
      
      shap_data_subsets <- list(list(split_from_data_frame$data, "full"), list(train, "train"), list(test, "test"))
      
      for (i in seq_along(shap_data_subsets)) {
        
        ## prep the data using the recipe
        shap_data <- recipes::prep(prepped_recipe, shap_data_subsets[[i]][[1]], fresh = TRUE) %>%
          recipes::juice() %>%
          dplyr::select(-feature_of_interest, -subject_id)
        assign(paste0("shap_data_", shap_data_subsets[[i]][[2]]), shap_data, envir = shap_plot_env)
        
        ## fit workflow to new datasets
        shap_model_object <- parsnip::fit(final_workflow, shap_data_subsets[[i]][[1]])
        ## extract the model itself
        shap_model_object <- workflows::extract_fit_parsnip(shap_model_object)
        
        ## shap safety checks!
        n_rows <- nrow(shap_data)
        n_cols <- ncol(shap_data)
        # try and calc a conservative nsim, otherwise choose 10.
        safe_nsim <- max(10, floor(1200000 / (n_rows * n_cols))) # 20 features x 300 samples x 200 sims = 1200000
        ## for smaller datasets, the above could lead to huge nsim. lets set max at 200.
        safe_nsim <- ifelse(safe_nsim > 199, 200, safe_nsim)
        
        ## warning if shap analysis looks like its going to take a long time
        if ((n_cols * n_rows) > 500000) {
          warning("This input dataset is pretty large for a SHAP analysis. This may take a long time, potentially exceeding walltime limits for shared resources (e.g., HPCs)", immediate. = T)
        }
        
        message(glue::glue("Running SHAP with nsim = {safe_nsim}"))
        
        ## start a parallel process
        cl <- parallel::makeForkCluster(as.numeric(parallel_workers))
        doParallel::registerDoParallel(cl)
        
        # Compute SHAP values
        shap_explanations <- fastshap::explain(
          object = shap_model_object,
          X = shap_data,
          pred_wrapper = pfun,
          nsim = safe_nsim,
          adjust = TRUE,
          parallel = TRUE
        )
        
        parallel::stopCluster(cl)
        
        assign(paste0("shap_explanations_", shap_data_subsets[[i]][[2]]), shap_explanations, envir = shap_plot_env)
        
        # SHAP object for plotting
        sv <- shapviz::shapviz(shap_explanations, X = shap_data)
        assign(paste0("sv_", shap_data_subsets[[i]][[2]]), sv, envir = shap_plot_env)
        
        ## append scores to logger
        log_shap_analysis(sv_object = sv, top_features = 10, data_subset = shap_data_subsets[[i]][[2]])
        
        # Generate and save plot
        plot <- shap_plot(
          sv = sv,
          label = label,
          data_subset_label = shap_data_subsets[[i]][[2]],
          split_from_data_frame = split_from_data_frame,
          filename = filename,
          output_dir = output_dir,
          data_subset_index = i,
          feature_type = feature_type
        )
        assign(paste0("plot_", shap_data_subsets[[i]][[2]]), plot, envir = shap_plot_env)
      }
      
      
    }, error = function(e) {
      shap.error.occured <<- TRUE
      error_message <<- e$message
      NULL
    })
  }
  
  # --- Save and return results ---
  if (shap.error.occured) {
    message(paste("SHAP analysis encountered an issue and all output files may not have been generated:", error_message))
    if (!is.null(error_message)) { 
      ## attempt to still return what was written to shap_plot_env
      save(list = ls(envir = shap_plot_env), 
           envir = shap_plot_env,
           file = file.path(paste0(output_dir, "/shap_inputs_", filename, ".RData")),
           compress = "gzip"
      )
      ## return error message
      message("Error: ", error_message)
    }
  } else {
    message("✅ SHAP analysis completed successfully.")
    save(list = ls(envir = shap_plot_env), 
         envir = shap_plot_env,
         file = file.path(paste0(output_dir, "/shap_inputs_", filename, ".RData")),
         compress = "gzip"
    )
  }
  
  # --- Clean up large local objects ---
  rm(list = ls(envir = shap_plot_env), envir = shap_plot_env)
  gc(verbose = FALSE)
  
  return(invisible(list(
    success = !shap.error.occured,
    shap_plot_env = shap_plot_env,
    error_message = error_message
  )))
  
}

shap_plot <- function(
    sv,
    label,
    data_subset_label,
    split_from_data_frame,
    filename,
    output_dir,
    data_subset_index,
    feature_type
) {
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Determine class labels (assumes binary classification)
  class_levels <- levels(as.factor(split_from_data_frame$data$feature_of_interest))
  if (length(class_levels) < 2) {
    stop("Insufficient factor levels for feature_of_interest.")
  }
  ## log shap viz positive class
  if (feature_type == "factor") {
    logger::log_info(paste0("Positive class (shap_vis): ", class_levels[1]))
  }
  
  # Create plot
  ## MODIFY THESE PARAMETERS IF YOU WANT THE PLOT TO LOOK DIFFERENTLY!!
  plot <- shapviz::sv_importance(
    sv,
    kind = "bee",
    show_numbers = TRUE,
    bee_width = 0.2,
    max_display = 10
  ) +
    ggtitle(label = paste0("SHAP: ", label, " (", data_subset_label, ")")) +
    labs(x = ifelse(feature_type == "factor", paste0(
      "predictive of ",
      class_levels[2],
      " < SHAP > predictive of ",
      class_levels[1]
    ), paste0(
      "low response < SHAP > high response "
    ))) +
    theme_bw(base_size = 14)
  
  # Construct filename and save plot
  filename_out <- file.path(output_dir, paste0("shap_", filename, "_", data_subset_label, ".pdf"))
  
  ggplot2::ggsave(
    plot = plot,
    filename = filename_out,
    width = pmax(0.1 * max(nchar(colnames(sv$X))), 6),
    height = 4.5,
    units = "in"
  )
  
  message("SHAP plot saved to: ", filename_out)
  
  return(plot)
}

## This is a helper script to shorten the long names of shap plots. It doesnt
## get used in this codebase, but its too good not to exist somewhere. 
## TODO: If we organize code, this should live with the shap_analysis() code
shap_shorten_colnames <- function(shap_sv_obj, splits) {
  # Function to extract the last matching split and everything after it
  shorten_name <- function(name, splits) {
    matches <- sapply(splits, function(split) {
      regexpr(split, name, fixed = TRUE)
    })
    
    # Filter valid matches (not -1), and get the last one
    valid_matches <- which(matches != -1)
    if (length(valid_matches) == 0) {
      return(name)  # no split found, return original
    }
    
    last_match_pos <- max(matches[valid_matches])
    substring(name, last_match_pos)
  }
  
  # Apply shortening to all column names
  new_colnames <- sapply(colnames(shap_sv_obj$X), shorten_name, splits = splits, USE.NAMES = FALSE)
  
  # Create new object with shortened column names
  obj_new <- shap_sv_obj
  colnames(obj_new$X) <- new_colnames
  colnames(obj_new$S) <- new_colnames
  
  return(obj_new)
}

log_shap_analysis <- function(sv_object, top_features, data_subset) {
  
  ## Add shap header to log
  logger::log_info(paste0("Mean absolute SHAP results: ", data_subset, " data (*up to top 10 features)"))
  
  ## get the top n features based on mean abs shap values
  tmp <- data.frame(feature=colnames(sv_object$S), 
                    mean_abs_shap=colMeans(abs(sv_object$S))) %>% 
    tibble::remove_rownames() %>% 
    dplyr::arrange(., desc(mean_abs_shap)) %>%
    dplyr::slice_head(., n = top_features)
  
  ## add them individually to logger
  for (n in seq(1:nrow(tmp))) {
    logger::log_info(paste0(tmp[n,1], " | ", tmp[n,2]))
  }
}

