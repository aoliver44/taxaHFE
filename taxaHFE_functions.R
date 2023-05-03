## HFE FUNCTIONS
## v0.3.0a.2

## read in microbiome data =====================================================
read_in_microbiome <- function(input) {
  ## read in txt, tsv, or csv microbiome data
  if (strsplit(basename(input), split="\\.")[[1]][2] %in% c("tsv","txt")) {
    suppressMessages(readr::read_delim(file = input, delim = "\t", skip = 0, name_repair = "minimal") %>% dplyr::select(., -any_of(c("NCBI_tax_id", "clade_taxid"))))
  } else {
    suppressMessages(readr::read_delim(file = input, delim = ",", skip = 0, name_repair = "minimal") %>% dplyr::select(., -any_of(c("NCBI_tax_id", "clade_taxid"))))
  }

}

## read in metadata  ===========================================================
read_in_metadata <- function(input, subject_identifier, label) {
  if (strsplit(basename(input), split="\\.")[[1]][2] %in% c("tsv","txt")) {
    suppressMessages(readr::read_delim(file = input, delim = "\t")) %>% 
      dplyr::select(., subject_identifier, label) %>% 
      dplyr::rename(., "subject_id" = subject_identifier) %>%
      rename(., "feature_of_interest" = label) %>%
      tidyr::drop_na()
  } else {
    suppressMessages(readr::read_delim(file = input, delim = ",")) %>% 
      dplyr::select(., subject_identifier, label) %>% 
      dplyr::rename(., "subject_id" = subject_identifier) %>%
      rename(., "feature_of_interest" = label) %>%
      tidyr::drop_na()
  }
}


## run safety checks!  =========================================================

mid_safety_checks <- function(type = opt$feature_type, label = opt$label, out = opt$output, input = hData, meta = metadata) {

  ## check if type was mis-specified
  if (type == "factor") {
    if(length(levels(as.factor(metadata$feature_of_interest))) > 9)
      stop("You are trying to predict 10 or more classes.\nThat is a bit much. Did you mean to do regression? (i.e., --feature_type numeric)")
  }
  
  ## check and see if output directory exists
  if (!dir.exists(dirname(out))) {
    stop("No output path found. Did you create the output directory?")
  }
  
  ## check and make sure colnames and metadata match
  input <- input %>% tibble::column_to_rownames(., var = "clade_name") %>% t() %>% as.data.frame()
  tmp_merge <- merge(meta, input, by.x = "subject_id", by.y = "row.names")
  
  if (NROW(tmp_merge) < 2) {
    stop("Your subject identifer doesnt match between input data and metadata.\n Or it's in the incorrect format.")
  }

}
## read in covariates  =========================================================

read_in_covariates <- function(input, subject_identifier) {
  if (strsplit(basename(input), split="\\.")[[1]][2] %in% c("tsv","txt")) {
    suppressMessages(readr::read_delim(file = input, delim = "\t")) %>% 
      dplyr::rename(., "subject_id" = subject_identifier) %>%
      tidyr::drop_na()
  } else {
    suppressMessages(readr::read_delim(file = input, delim = ",")) %>% 
      dplyr::rename(., "subject_id" = subject_identifier) %>%
      tidyr::drop_na()
  }
}

## convert to metaphlan  =======================================================

convert_to_hData <- function(input) {
  
  ## long vector of possible levels in hierarchical data
  levels <- c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "L11", "L12", "L13", "L14", "L15")
  
  ## count number of splits (+1) in hierarchical data by "|" symbol
  num_levels <- max(stringr::str_count(input$clade_name, "\\|"))
  levels <- levels[1:(num_levels + 1)]
  assign(x = "levels", value = levels, envir = .GlobalEnv)
  
  cat("\n\n", "Attempting to convert to metaphlan-style input...")
  
  ## split raw data by pipe symbol into number of expected parts
  input <- input %>% 
    dplyr::relocate(., clade_name) %>%
    tidyr::separate(., col = clade_name, into = levels, sep = "\\|", extra = "merge")

  ## summarize each level into its own dataframe
  ## use ultra fast data.table processesing
  for (level in seq(length(levels))) {
    
    summarize_level <- paste0("hData_L", level)
    hData_tmp <- input %>% 
      dplyr::select(., 1:level, (length(levels) + 1):dplyr::last_col()) %>% 
      tidyr::unite(., "clade_name", L1:levels[level], sep = "|", remove = F, na.rm = T) %>%
      dplyr::select(., -c( L1:levels[level])) %>%
      dtplyr::lazy_dt() 
    hData_tmp <- hData_tmp %>% 
      group_by(., clade_name) %>% 
      summarise_all(base::sum) %>% 
      as.data.frame() %>% 
      dplyr::filter(., !grepl(pattern = "NA", clade_name))
    assign(summarize_level, hData_tmp, envir = .GlobalEnv)
    
  }
}


## apply filters  ==============================================================

apply_filters <- function(input) {
  
  ## this CV check checks to make sure everything is going to get roughly
  ## divided by the same number for the abundance filter
  ## if not, then it will break.
  CV <- function(x){
    (sd(x)/mean(x))*100
  }
  
  cv_check <- input %>% 
    dplyr::filter(., !grepl("\\|", clade_name)) %>% 
    tibble::column_to_rownames(., var = "clade_name") %>% 
    dplyr::summarise_all(sum) %>% t() %>% as.data.frame() %>% dplyr::pull()
  
  if ((CV(cv_check) >= 0.1) == TRUE) {
    cat("\n\n\n", "###########################\n", "ERROR: CV too large \n", "###########################")
    cat("\n", "Program is stopping because your data is not normalized. It should be rarefied or in relative abundance.")
    stop()
  }
  
  cat("\n\n", "###########################\n", "Applying filtering steps...\n", "###########################")
  post_metaphlan_transformation_feature_count <- NROW(input)
  ## remove taxa/rows that are 99% zeros (1% prevalence filter)
  input <- input[rowSums(input[,2:NCOL(input)] == 0) <= (NCOL(input[,2:NCOL(input)])*0.99), ]
  cat("\n\n Prevelance filter: ")
  prev_filter <- NROW(input)
  assign(x = "prev_filter", value = prev_filter, envir = .GlobalEnv)
  cat(paste0(round((((post_metaphlan_transformation_feature_count - prev_filter)/ post_metaphlan_transformation_feature_count ) * 100)), "% of features dropped due to 1% prevelance filter.\n"))
  
  ## Remove very low abundant features ===========================================
  ## remove taxa/rows that are below 0.0001 relative abundance
  cat("\n Low abundance filter: ")
  hData_mean_total_abundance <- input %>% 
    dplyr::filter(., !grepl("\\|", clade_name)) %>% 
    tibble::column_to_rownames(., var = "clade_name") %>% 
    summarise_all(sum) %>% rowMeans()
  
  
  ## abundance filter
  hData_abund_filter <- input %>% tibble::remove_rownames() %>% tibble::column_to_rownames(., var = "clade_name")
  hData_abund_filter$mean_abundance <- rowMeans(hData_abund_filter)
  hData_abund_filter <- hData_abund_filter %>% 
    dplyr::filter(., (mean_abundance / hData_mean_total_abundance) >= 0.0001) %>% 
    tibble::rownames_to_column(., var = "clade_name") %>% 
    dplyr::pull(., clade_name)
  
  hData <- input %>% dplyr::filter(., clade_name %in% hData_abund_filter)
  assign(x = "hData", value = hData, envir = .GlobalEnv)
  cat(paste0(round((((prev_filter - NROW(hData)) / prev_filter) * 100)), "% of post-prevalence filtered features dropped due \n to abundance filter (rel. abund. >10e-4) \n\n"))
  cat(NROW(hData), "taxa retained for downstream analysis...\n")
}

## make make_taxa_split dataframe ============================================== 
make_taxa_split_df <- function(input) {
  
  ## long vector of possible levels in hierarchical data
  levels <- c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "L11", "L12", "L13", "L14", "L15")
  
  ## count number of splits (+1) in hierarchical data by "|" symbol
  num_levels <- max(stringr::str_count(input$clade_name, "\\|"))
  levels <- levels[1:(num_levels + 1)]
  
  ## split raw data by pipe symbol into number of expected parts
  taxa_only_split <- input %>% 
    dplyr::relocate(., clade_name) %>%
    tidyr::separate(., col = clade_name, into = levels, sep = "\\|", extra = "merge", remove = FALSE) %>%
    dplyr::select(., 1:(length(levels) + 1)) 
  
  ## add some columns to this dataframe (full clade name, full taxa abundance, number of NAs)
  ## NOTE: future scripts could add an abundance threshold filter here
  taxa_only_split$taxa_abundance <- rowSums(hData[,2:NCOL(hData)])
  taxa_only_split$na_count <- rowSums(is.na(taxa_only_split))
  
  assign("taxa_only_split", taxa_only_split, envir = .GlobalEnv)
  assign("levels", levels, envir = .GlobalEnv)

}

## write separate files to test summary levels =================================

write_summary_files <- function(input, output) {
  
  ## long vector of possible levels in hierarchical data
  levels <- c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "L11", "L12", "L13", "L14", "L15")
  
  ## count number of splits (+1) in hierarchical data by "|" symbol
  num_levels <- max(stringr::str_count(input$clade_name, "\\|"))
  levels <- levels[1:(num_levels + 1)]
  
  ## split raw data by pipe symbol into number of expected parts
  hData_summary <- input %>% 
    dplyr::relocate(., clade_name) %>%
    tidyr::separate(., col = clade_name, into = levels, sep = "\\|", extra = "merge", remove = FALSE) %>%
    dplyr::mutate(., level = ((stringr::str_count(hData$clade_name, "\\|")) + 1))
  
  count = 1
  for (i in seq(length(levels))) {
    
    ## select different levels and write them to file
    file_summary <- hData_summary %>% 
      dplyr::filter(., level == i) %>%
      dplyr::select(., -level, -dplyr::any_of(levels)) %>%
      tibble::column_to_rownames(., var = "clade_name") %>%
      t() %>%
      as.data.frame() %>%
      tibble::rownames_to_column(., var = "subject_id") %>%
      janitor::clean_names()

    file_merge <- merge(metadata, file_summary, by = "subject_id")
    
    filename <- paste0("_level_", count, ".csv")
    readr::write_csv(x = file_merge, file = paste0(tools::file_path_sans_ext(output), filename))
    count = count + 1
  }
  
}

## function to calculate class imbalance =======================================

calc_class_frequencies <- function(input = metadata, feature_type = opt$feature_type) {
  
  if (feature_type == "factor") {
    ## create a vector to deal with class imbalance
    ## this is infomred by: https://github.com/imbs-hl/ranger/issues/167
    ## note in the ranger model, they use sample.fraction and replace = T
    class_frequencies <- metadata %>% 
      dplyr::count(feature_of_interest) %>% 
      dplyr::mutate(prop = prop.table(n)) %>% dplyr::pull(prop)
    ## make the class frequencies a fraction of the entire data to help
    ## prevent overfitting. Janky, but i think this is better than nothing
    ## edit...im gonna leave this in here, but since ranger uses OOB to calc
    ## scores maybe its best to use all the data. Not ideal, but because
    ## sample sizes are usually a little small i think its better to use all
    ## data??
    class_frequencies <- class_frequencies * as.numeric(opt$subsample)
    
    assign("class_frequencies", class_frequencies, envir = .GlobalEnv)
  }
  
}

## main competition function ===================================================

taxaHFE_competition <- function(input = hData, feature_type = opt$feature_type, cores = 4, output) {
  
  cat("\n\n", "##################################\n", "Starting hierarchical competitions\n", "##################################\n\n")
  
  ## helper function
  `%!in%` <- Negate(`%in%`)
  
  ## lists to loop over for parents and children levels
  parents <- purrr::discard(dplyr::lead(rev(levels)), is.na)
  children <- rev(levels)
  count = 1
  
  ## correlation vs RF list
  cor_vs_rf <- c()
  
  for (parent in parents) {
    
    ## make a dataframe of all the species and how many sub-species each species has
    children_counts <- taxa_only_split %>%
      dplyr::group_by(., .data[[parent]]) %>%
      dplyr::summarize(., count = n_distinct(.data[[children[count]]])) %>%
      tidyr::drop_na()
    
    ## start progress bar
    pb <- progress_bar$new( format = " Child vs Parents [:bar] :percent in :elapsed", total = length(children_counts[children_counts$count > 1, ][,1] %>% dplyr::pull()), clear = FALSE, width= 60) 
    
    ## start a counter to keep track of each for-loop 
    sub_count = 1
    
    ## loop through all the species with greater than 1 sub-species. If species
    ## has only itself, it will automatically be kept and used in the genus step.
    for (parent_trial in children_counts[children_counts$count > 1, ][,1] %>% dplyr::pull()) {
      assign("parent_trial", parent_trial, envir = .GlobalEnv)
      ## progress bar tick
      pb$tick()
      ## create a dataframe of the parent and children. In this case the abundance
      ## of a species and the subspecies across samples.
      hData_parent <- input %>%
        dplyr::filter(., grepl(pattern = paste0("\\b", parent_trial, "\\b"), clade_name)) %>%
        tidyr::separate(., col = clade_name, into = levels, sep = "\\|") %>%
        dplyr::select(., .data[[children[count]]], where(is.numeric)) %>%
        dplyr::rename(., "child" = 1) %>% 
        tidyr::replace_na(., list(child = "PARENT")) %>%
        dplyr::group_by(., child) %>%
        dplyr::summarise(., across(where(is.numeric), ~ sum(.))) %>% 
        tibble::column_to_rownames(., var = "child") %>%
        t() %>%
        as.data.frame() 
      
      ## if no columns named parent exists (species renamed parent in previous step)
      ## create column by summing up children (eg,sub-species)
      if ("PARENT" %!in% colnames(hData_parent)) { 
        hData_parent$PARENT <- rowSums(hData_parent)
      }
      ## merge with metadata
      hData_parent_merge <- merge(metadata, hData_parent, by.x = "subject_id", by.y = "row.names")
      
      ### CORRELATION #####
      
      ## correlate parent with children. If parent (eg, species) is highly (pearson = 0.95) 
      ## correlated with child, drop the highly correlated child...they dont bring
      ## more information to the table that is otherwise carried in the parent (species in this case)
      cor_drop <- suppressMessages(corrr::correlate(hData_parent)) %>% 
        corrr::focus(., PARENT) %>% dplyr::filter(., PARENT >= as.numeric(opt$cor_level)) %>% 
        dplyr::pull(., term)
      
      hData_parent_merge_cor_subset <- hData_parent_merge %>% 
        dplyr::select(., -all_of(cor_drop)) %>%
        select(., -subject_id)
      
      ## if, after dropping highly correlated children, only parent is left,
      ## drop from taxa_only_split (dataframe keeping track of all kept features)
      ## and move on to the next species...ELSE, run a random forest, below
      if (NCOL(hData_parent_merge_cor_subset) < 3) {
        ## so this is selecting the features in taxa_only_split where the for-loop
        ## variable is found in the Parent column (in this case species), and there are no
        ## NAs (which is expected, because the hierarchy is complete if were looking at
        ## species vs subspecies. If we were looking at genus vs species, we'd expect 1 NA, because
        ## after this step, only a species (PARENT) or a subspecies(CHILD) moved forward, and one was 
        ## dropped...thus a NA was created. Just a way to make the dropping more selective.)
        mutate_child = children[count]
        taxa_only_split <- taxa_only_split %>% dplyr::mutate(., !!mutate_child := ifelse((.data[[parents[count]]] == parent_trial & !is.na(.data[[parents[count]]]) & na_count == (count - 1)), "drop_dis", .data[[children[count]]]))
        taxa_only_split <- taxa_only_split %>% dplyr::filter(., !grepl(pattern = "drop_dis", x = .data[[children[count]]]))
        assign("taxa_only_split", taxa_only_split, envir = .GlobalEnv)
        ## So some subspecies (CHILDREN) were not correlated with the species (PARENT)...
        ## Do these species bring more information to the table with regards to the
        ## feature_of_interest?
        
        ## keep tally of how many features were collapsed due to correlation vs RF
        cor_vs_rf <- append(cor_vs_rf, 1, after = length(cor_vs_rf))
        
      } else {
        
        ## RF MODEL #####
        
        ## keep tally of how many features were collapsed due to correlation vs RF
        cor_vs_rf <- append(cor_vs_rf, 0, after = length(cor_vs_rf))
        
        ## if the feature of interest is a factor, do a RF Classification
        if (feature_type == "factor") {
          
          ## create an intial model and keep track of the variable.importance
          ## this will help us decide if the PARENT (species) or the CHILD (sub-speceies)
          ## brings more information to the table with regards to feature_of_interest
          model <- ranger::ranger(as.factor(feature_of_interest) ~ . , data = hData_parent_merge_cor_subset, importance = "impurity_corrected", seed = 42, sample.fraction = class_frequencies, replace = TRUE, num.threads = as.numeric(cores))
          model_importance <- as.data.frame(model$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
          
          ## welp, the initial model miiight be correct, but lets permute that process
          ## by looping over some random seeds (10, as set in a very early variable nperm at the top)
          ## and keeping track of the variable.importance. We can then average that and have
          ## a more sure guess whether a Parent or Child is more important with regards to 
          ## the feature_of_interest
          for (seed in sample(1:1000, nperm)) {
            model_tmp <- ranger::ranger(as.factor(feature_of_interest) ~ . , data = hData_parent_merge_cor_subset, importance = "impurity_corrected", seed = seed, sample.fraction = class_frequencies, replace = TRUE, num.threads = as.numeric(cores))
            model_importance_tmp <- as.data.frame(model_tmp$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
            suppressWarnings(model_importance <- merge(model_importance, model_importance_tmp, by = "taxa"))
            
          }
          colnames(model_importance)[2:(nperm + 2)] <- paste0("permutation_", seq(1,nperm + 1))
          model_importance$average <- rowMeans(model_importance[, 2:(nperm + 2)])
          
          ## this else statement does the sample as the above few lines, just for a continous
          ## feature_of_interest...with RF Regression.
        } else {
          model <- ranger::ranger(as.numeric(feature_of_interest) ~ . , data = hData_parent_merge_cor_subset, importance = "impurity_corrected", seed = 42, sample.fraction = 0.7, num.threads = as.numeric(cores))
          model_importance <- as.data.frame(model$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
          for (seed in sample(1:1000, nperm)) {
            model_tmp <- ranger::ranger(as.numeric(feature_of_interest) ~ . , data = hData_parent_merge_cor_subset, importance = "impurity_corrected", seed = seed, sample.fraction = 0.7, num.threads = as.numeric(cores))
            model_importance_tmp <- as.data.frame(model_tmp$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
            suppressWarnings(model_importance <- merge(model_importance, model_importance_tmp, by = "taxa"))
            
          }
          colnames(model_importance)[2:(nperm + 2)] <- paste0("permutation_", seq(1,nperm + 1))
          model_importance$average <- rowMeans(model_importance[, 2:(nperm + 2)])
          
          
        }
        
        ### RF WINNER - PARENT ####
        ## Ok we have the permuted variable.importance. Do, via the average, any
        ## children (sub-species) beat the parents with regards to discriminatory information
        ## with regards to the feature_of_interest? 
        ## If the largest average variable importance is PARENT:
        if ((model_importance %>% dplyr::arrange(., desc(average)) %>% slice_head(., n = 1) %>% dplyr::pull(., taxa)) == "PARENT") {
          ## drop children from master taxa file taxa_only_split
          mutate_child = children[count]
          taxa_only_split <- taxa_only_split %>% dplyr::mutate(., !!mutate_child := ifelse((.data[[parents[count]]] == parent_trial & !is.na(.data[[parents[count]]]) & na_count == (count - 1)), "drop_dis", .data[[children[count]]]))
          taxa_only_split <- taxa_only_split %>% dplyr::filter(., !grepl(pattern = "drop_dis", x = .data[[children[count]]]))
          assign("taxa_only_split", taxa_only_split, envir = .GlobalEnv)
          
        } else { 
          ### RF WINNER - CHILD ####
          mutate_child = children[count]
          parent_importance <- model_importance$average[model_importance$taxa == "PARENT"]
          
          children_toss <- model_importance %>% dplyr::filter(., average < parent_importance) %>% dplyr::pull(., taxa)
          children_toss_zero <- model_importance%>% dplyr::filter(., average < 0) %>% dplyr::pull(., taxa) 
          children_toss <- unique(c(children_toss, cor_drop, children_toss_zero))
          
          ## drop parent
          mutate_at_list <- rev(parents)[1:(length(parents) + 1 - count)]
          taxa_only_split <- taxa_only_split %>% 
            dplyr::mutate(., !!mutate_child := ifelse((.data[[parents[count]]] == parent_trial & !is.na(.data[[parents[count]]]) & is.na(.data[[children[count]]])), "drop_dis", .data[[children[count]]]))
          taxa_only_split <- taxa_only_split %>% 
            mutate_at(mutate_at_list, ~ 
                        replace(., .data[[parents[count]]] == parent_trial, NA)) %>%
            dplyr::filter(., !grepl(pattern = "drop_dis", x = .data[[children[count]]]))
          
          ## drop children that didnt win against parent
          taxa_only_split <- taxa_only_split %>% 
            dplyr::mutate(., !!mutate_child := ifelse((.data[[children[count]]] %in% children_toss & is.na(.data[[parents[count]]]) & !is.na(.data[[children[count]]])), "drop_dis", .data[[children[count]]]))
          taxa_only_split <- taxa_only_split %>% 
            dplyr::filter(., !grepl(pattern = "drop_dis", x = .data[[children[count]]]))
          assign("taxa_only_split", taxa_only_split, envir = .GlobalEnv)
          
          
        }
      }
      ### PROGRESS ####
      if (sub_count == length(children_counts[children_counts$count > 1, ][,1] %>% dplyr::pull())) message(paste0("  ", "Done with ", parents[count], " versus ", children[count], ". Moving on to next levels..."))
      sub_count = sub_count + 1
    }
    
    hData <- hData %>%
      dplyr::filter(., clade_name %in% taxa_only_split$clade_name)
    assign("hData", hData, envir = .GlobalEnv)
    count = count + 1
    
  }
  
  ## write no SF files
  hData_output <- hData %>%
    tibble::column_to_rownames(., var = "clade_name") %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(., var = "subject_id") %>%
    janitor::clean_names()
  output_nsf <- merge(metadata, hData_output, by.x = "subject_id", by.y = "subject_id")
  assign("output_nsf", output_nsf, envir = .GlobalEnv)
  assign("output_nsf_count", (NCOL(output_nsf) - 2), envir = .GlobalEnv)
  readr::write_delim(file = paste0(tools::file_path_sans_ext(output), "no_sf.txt"), x = output_nsf, delim = "\t")
  readr::write_delim(file = paste0(tools::file_path_sans_ext(output), "_taxa_list_no_sf.txt"), x = taxa_only_split, delim = "\t")
  
  ## assign cor vs rf battle tally to global env
  assign("cor_vs_rf", cor_vs_rf, envir = .GlobalEnv)
  
  ## print results of cor vs rf battle tally
  cat(paste0("\nPercent of parent-child competitions decided through correlation: ", round((sum(cor_vs_rf)/ length(cor_vs_rf) * 100), digits = 2), "\n"))
  cat(paste0("Percent of parent-child competitions decided through RF models: ", round(((length(cor_vs_rf) - sum(cor_vs_rf))/ length(cor_vs_rf) * 100), digits = 2), "\n"))
}

## super filter ================================================================

super_filter <- function(input = hData, feature_type = "factor", cores = 4, subsample = opt$subsample, output) {
  
  cat("\n\n", "##################################\n", "Starting Super-filter\n", "##################################\n\n")
  
  hData_sf <- input %>%
    tibble::column_to_rownames(., var = "clade_name") %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(., var = "subject_id") %>%
    janitor::clean_names()
  
  hData_sf <- merge(metadata, hData_sf, by.x = "subject_id", by.y = "subject_id")
  output_sf <- hData_sf %>% dplyr::select(., -subject_id) 
  nperm = nperm + 190
  
  pb <- progress_bar$new( format = " Super-filter [:bar] :percent in :elapsed", total = nperm, clear = FALSE, width= 60)
  
  if (feature_type == "factor") {  
    model <- ranger::ranger(as.factor(feature_of_interest) ~ . , data = output_sf, importance = "impurity_corrected", seed = 42, sample.fraction = class_frequencies, replace = TRUE, num.threads = as.numeric(cores))
    model_importance <- as.data.frame(model$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
    for (seed in sample(1:1000, nperm)) {
      
      model_tmp <- ranger::ranger(as.factor(feature_of_interest) ~ . , data = output_sf, importance = "impurity_corrected", seed = seed, sample.fraction = class_frequencies, replace = TRUE, num.threads = as.numeric(cores))
      model_importance_tmp <- as.data.frame(model_tmp$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
      suppressWarnings(model_importance <- merge(model_importance, model_importance_tmp, by = "taxa"))
      
      pb$tick()
    }
    colnames(model_importance)[2:(nperm + 2)] <- paste0("permutation_", seq(1,nperm + 1))
    model_importance$average <- rowMeans(model_importance[, 2:(nperm + 2)])
    model_importance <- model_importance %>% dplyr::relocate(., average)
    assign("model_importance", model_importance, envir = .GlobalEnv)
    
  } else {
    model <- ranger::ranger(as.numeric(feature_of_interest) ~ . , data = output_sf, importance = "impurity_corrected", seed = 42, sample.fraction = as.numeric(subsample), num.threads = as.numeric(cores))
    model_importance <- as.data.frame(model$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
    for (seed in sample(1:1000, nperm)) {
      
      model_tmp <- ranger::ranger(as.numeric(feature_of_interest) ~ . , data = output_sf, importance = "impurity_corrected", seed = seed, sample.fraction = as.numeric(subsample), num.threads = as.numeric(cores))
      model_importance_tmp <- as.data.frame(model_tmp$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
      suppressWarnings(model_importance <- merge(model_importance, model_importance_tmp, by = "taxa"))
      
      pb$tick()
    }
    colnames(model_importance)[2:(nperm + 2)] <- paste0("permutation_", seq(1,nperm + 1))
    model_importance$average <- rowMeans(model_importance[, 2:(nperm + 2)])
    model_importance <- model_importance %>% dplyr::relocate(., average)
    assign("model_importance", model_importance, envir = .GlobalEnv)
  }
  
  
  model_importance_list <- model_importance %>% dplyr::filter(., average > mean(average)) %>% dplyr::filter(., average > 0) %>% dplyr::pull(., taxa)
  output_sf <- hData_sf %>% dplyr::select(., subject_id, feature_of_interest, all_of(model_importance_list))
  taxa_only_split$clade_name <- janitor::make_clean_names(taxa_only_split$clade_name)
  taxa_only_split <- taxa_only_split %>% dplyr::filter(., clade_name %in% model_importance_list)
  
  assign("output_sf", output_sf, envir = .GlobalEnv)
  assign("taxa_only_split_sf", taxa_only_split, envir = .GlobalEnv)
  
  readr::write_delim(file = output, x = output_sf, delim = "\t")
  readr::write_delim(file = paste0(tools::file_path_sans_ext(output), "_taxa_list.txt"), x = taxa_only_split, delim = "\t")
  
  cat("\n\n Features (no super filter): ", output_nsf_count)
  cat("\n Features (super filter): ", (NCOL(output_sf) - 2), "\n\n")
  
}

## write figures ===============================================================

write_figure <- function(input, output) {
  
  top_features <- model_importance %>% dplyr::filter(., average > mean(average)) %>% 
    dplyr::filter(., average > 0) %>% dplyr::arrange(., desc(average)) %>% dplyr::pull(., taxa)
  
  top_features <- top_features[1:pmin(10, length(top_features))]
  
  figure_data <- output_sf %>% dplyr::select(., feature_of_interest, any_of(top_features)) %>%
    reshape2::melt(id.vars = "feature_of_interest")
  
  if (opt$feature_type == "factor") {
    suppressWarnings(ggplot(data = figure_data) +
                       aes(x = as.factor(feature_of_interest), y = log(value)) +
                       geom_boxplot(aes(fill = as.factor(feature_of_interest)), outlier.alpha = 0) +
                       geom_point(position = position_jitter(width = 0.2), alpha = 0.4) +
                       facet_wrap( ~ variable, scales = "free_y", ncol = 1, labeller = labeller(groupwrap = label_wrap_gen(10))) +
                       theme_bw() +
                       theme(strip.text.x = element_text(size = 7), legend.position = "none", text = element_text(color = "black")) + 
                       labs(y = "ln(Relative abundance)", x = "Feature of Interest") +
                       ggsci::scale_fill_jama() + coord_flip())
    
    ggsave(filename = paste0(tools::file_path_sans_ext(output), "_plot.pdf"), device = "pdf", dpi = "retina", width = pmax((max(nchar(as.character(figure_data$variable))) * 0.0555), 3), height = 10, units = "in")
    
  } else {
    suppressWarnings(ggplot(data = figure_data) +
                       aes(x = feature_of_interest, y = log(value)) +
                       geom_point() +
                       geom_smooth(method = "lm") +
                       facet_wrap( ~ variable, scales = "free_y") +
                       theme_bw() + theme(strip.text.x = element_text(size = 5)))
    
    ggsave(filename = paste0(tools::file_path_sans_ext(output), "_plot.pdf"), device = "pdf", dpi = "retina", width = 12, height = 8, units = "in")
    
  }
  
  
}

## write files for old_HFE =====================================================

write_old_hfe <- function(input = hData, output) {
  
  ## long vector of possible levels in hierarchical data
  levels <- c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "L11", "L12", "L13", "L14", "L15")
  
  ## count number of splits (+1) in hierarchical data by "|" symbol
  num_levels <- max(stringr::str_count(input$clade_name, "\\|"))
  levels <- levels[1:(num_levels + 1)]

  
  input <- input %>% 
    dplyr::relocate(., clade_name) %>%
    tidyr::separate(., col = clade_name, into = levels, sep = "\\|", extra = "merge")
  
  

  ## split raw data by pipe symbol into number of expected parts
  if ("L7" %in% colnames(input)) {
    taxonomy <- input %>% dplyr::select(., L1:L7)
    input <- input %>% 
      dplyr::select(., L7, where(is.numeric)) %>%
      dplyr::group_by(., L7) %>%
      dplyr::summarise_if(is.numeric, sum)
    taxonomy <- taxonomy[taxonomy$L7 %in% input$L7, ]
    taxonomy <- taxonomy %>% tidyr::drop_na() %>% dplyr::filter(., !duplicated(L7))
    input_taxa_merge <- merge(taxonomy, input, by = "L7")
    
    input_taxa_merge$index <- (1001:(NROW(input_taxa_merge) + 1000))
    input_taxa_merge$L1 <- "k__Bacteria"
    
    input_taxa_merge <- input_taxa_merge %>%
      dplyr::relocate(., index, L1,L2,L3,L4,L5,L6,L7)
    
    readr::write_delim(x = input_taxa_merge[1:8], file = paste0(tools::file_path_sans_ext(output), "_old_hfe_taxa.txt"), col_names = FALSE, delim = "\t")
    readr::write_delim(x = input_taxa_merge %>% dplyr::select(., 1,9:dplyr::last_col()), file = paste0(tools::file_path_sans_ext(output), "_old_hfe_otu.txt"), col_names = FALSE, delim = "\t")
    
    metadata_order <- colnames(input_taxa_merge[,9:NCOL(input_taxa_merge)])
    
    metadata_list <- metadata %>% dplyr::arrange(match(subject_id, metadata_order)) %>%
      pull(., feature_of_interest)
    
    metadata_list <- as.data.frame(c("label", metadata_list))
    readr::write_delim(x = as.data.frame(t(metadata_list)), file = paste0(tools::file_path_sans_ext(output), "_old_hfe_label.txt"), col_names = FALSE, delim = "\t")
    
    
  } else {
    taxonomy <- input %>% dplyr::select(., L1:L6)
    input <- input %>% 
      dplyr::select(., L6, where(is.numeric)) %>%
      dplyr::group_by(., L6) %>%
      dplyr::summarise_if(is.numeric, sum)
    taxonomy <- taxonomy[taxonomy$L6 %in% input$L6, ]
    taxonomy <- taxonomy %>% tidyr::drop_na() %>% dplyr::filter(., !duplicated(L6))
    input_taxa_merge <- merge(taxonomy, input, by = "L6")
    
    input_taxa_merge$index <- (1001:(NROW(input_taxa_merge) + 1000))
    input_taxa_merge$L1 <- "k__Bacteria"
    
    input_taxa_merge <- input_taxa_merge %>%
      dplyr::relocate(., index, L1,L2,L3,L4,L5,L6)
    
    readr::write_delim(x = input_taxa_merge[1:7], file = paste0(tools::file_path_sans_ext(output), "_old_hfe_taxa.txt"), col_names = FALSE, delim = "\t")
    readr::write_delim(x = input_taxa_merge %>% dplyr::select(., 1,8:dplyr::last_col()), file = paste0(tools::file_path_sans_ext(output), "_old_hfe_otu.txt"), col_names = FALSE, delim = "\t")
    
    metadata_order <- colnames(input_taxa_merge[,9:NCOL(input_taxa_merge)])
    
    metadata_list <- metadata %>% dplyr::arrange(match(subject_id, metadata_order)) %>%
      pull(., feature_of_interest)
    
    metadata_list <- as.data.frame(c("label", metadata_list))
    readr::write_delim(x = as.data.frame(t(metadata_list)), file = paste0(tools::file_path_sans_ext(output), "_old_hfe_label.txt"), col_names = FALSE, delim = "\t")
    
  }
}

## competition function with covariates ========================================
## essentially a copy of main competition function just with added
## merge and consideration of additional features

taxaHFE_competition_covariates <- function(input = hData, covariates, feature_type = opt$feature_type, cores = 4, output) {
  
  cat("\n\n", "####################################\n", "Starting competition with covariates\n", "####################################\n\n")
  
  ## helper function
  `%!in%` <- Negate(`%in%`)

  ## lists to loop over for parents and children levels
  parents <- purrr::discard(dplyr::lead(rev(levels)), is.na)
  children <- rev(levels)
  count = 1
  
  ## correlation vs RF list
  cor_vs_rf <- c()
  
  for (parent in parents) {
    
    ## make a dataframe of all the species and how many sub-species each species has
    children_counts <- taxa_only_split %>%
      dplyr::group_by(., .data[[parent]]) %>%
      dplyr::summarize(., count = n_distinct(.data[[children[count]]])) %>%
      tidyr::drop_na()
    
    ## start progress bar
    pb <- progress_bar$new( format = " Child vs Parents [:bar] :percent in :elapsed", total = length(children_counts[children_counts$count > 1, ][,1] %>% dplyr::pull()), clear = FALSE, width= 60) 
    
    ## start a counter to keep track of each for-loop 
    sub_count = 1
    
    ## loop through all the species with greater than 1 sub-species. If species
    ## has only itself, it will automatically be kept and used in the genus step.
    for (parent_trial in children_counts[children_counts$count > 1, ][,1] %>% dplyr::pull()) {
      assign("parent_trial", parent_trial, envir = .GlobalEnv)
      ## progress bar tick
      pb$tick()
      ## create a dataframe of the parent and children. In this case the abundance
      ## of a species and the subspecies across samples.
      hData_parent <- input %>%
        dplyr::filter(., grepl(pattern = paste0("\\b", parent_trial, "\\b"), clade_name)) %>%
        tidyr::separate(., col = clade_name, into = levels, sep = "\\|") %>%
        dplyr::select(., .data[[children[count]]], where(is.numeric)) %>%
        dplyr::rename(., "child" = 1) %>% 
        tidyr::replace_na(., list(child = "PARENT")) %>%
        dplyr::group_by(., child) %>%
        dplyr::summarise(., across(where(is.numeric), ~ sum(.))) %>% 
        tibble::column_to_rownames(., var = "child") %>%
        t() %>%
        as.data.frame() 
      
      ## if no columns named parent exists (species renamed parent in previous step)
      ## create column by summing up children (sub-species)
      if ("PARENT" %!in% colnames(hData_parent)) { 
        hData_parent$PARENT <- rowSums(hData_parent)
      }
      ## merge with metadata
      hData_parent_merge <- merge(metadata, hData_parent, by.x = "subject_id", by.y = "row.names")
    
      ### CORRELATION #####
      
      ## correlate parent with children. If parent (species) is highly (pearson = 0.95) 
      ## correlated with child, drop the highly correlated child...they dont bring
      ## more information to the table that is otherwise carried in the parent (species in this case)
      cor_drop <- suppressMessages(corrr::correlate(hData_parent)) %>% 
        corrr::focus(., PARENT) %>% dplyr::filter(., PARENT >= as.numeric(opt$cor_level)) %>% 
        dplyr::pull(., term)
      
      hData_parent_merge_cor_subset <- hData_parent_merge %>% 
        dplyr::select(., -all_of(cor_drop)) 

      ## if, after dropping highly correlated children, only species is left,
      ## drop from taxa_only_split (dataframe keeping track of all kept features)
      ## and move on to the next species...ELSE, run a random forest, below
      if (NCOL(hData_parent_merge_cor_subset) < 4) {
        ## so this is selecting the features in taxa_only_split where the for-loop
        ## variable is found in the Parent column (in this case species), and there are no
        ## NAs (which is expected, because the hierarchy is complete if were looking at
        ## species vs subspecies. If we were looking at genus vs species, we'd expect 1 NA, because
        ## after this step, only a species (PARENT) or a subspecies(CHILD) moved forward, and one was 
        ## dropped...thus a NA was created. Just a way to make the dropping more selective.)
        mutate_child = children[count]
        taxa_only_split <- taxa_only_split %>% dplyr::mutate(., !!mutate_child := ifelse((.data[[parents[count]]] == parent_trial & !is.na(.data[[parents[count]]]) & na_count == (count - 1)), "drop_dis", .data[[children[count]]]))
        taxa_only_split <- taxa_only_split %>% dplyr::filter(., !grepl(pattern = "drop_dis", x = .data[[children[count]]]))
        assign("taxa_only_split", taxa_only_split, envir = .GlobalEnv)
        ## So some subspecies (CHILDREN) were not correlated with the species (PARENT)...
        ## Do these species bring more information to the table with regards to the
        ## feature_of_interest?
        
        ## keep tally of how many features were collapsed due to correlation vs RF
        cor_vs_rf <- append(cor_vs_rf, 1, after = length(cor_vs_rf))
        
      } else {
        
        ## RF MODEL #####
        
        ## keep tally of how many features were collapsed due to correlation vs RF
        cor_vs_rf <- append(cor_vs_rf, 0, after = length(cor_vs_rf))
        
        ## merge with covariates
        hData_parent_merge_cor_subset <- merge(covariates, hData_parent_merge_cor_subset, by = "subject_id")
        hData_parent_merge_cor_subset <- hData_parent_merge_cor_subset %>% 
          dplyr::select(., -dplyr::any_of(opt$subject_identifier)) %>% 
          tidyr::drop_na()
        
        ## if the feature of interest is a factor, do a RF Classification
        if (feature_type == "factor") {
          ## create an intial model and keep track of the variable.importance
          ## this will help us decide if the PARENT (species) or the CHILD (sub-speceies)
          ## brings more information to the table with regards to feature_of_interest
          model <- ranger::ranger(as.factor(feature_of_interest) ~ . , data = hData_parent_merge_cor_subset, importance = "impurity_corrected", seed = 42, sample.fraction = 0.7, num.threads = as.numeric(cores))
          model_importance <- as.data.frame(model$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
          
          ## welp, the initial model miiight be correct, but lets permute that process
          ## by looping over some random seeds (10, as set in a very early variable nperm at the top)
          ## and keeping track of the variable.importance. We can then average that and have
          ## a more sure guess whether a Parent or Child is more important with regards to 
          ## the feature_of_interest
          for (seed in sample(1:1000, nperm)) {
            model_tmp <- ranger::ranger(as.factor(feature_of_interest) ~ . , data = hData_parent_merge_cor_subset, importance = "impurity_corrected", seed = seed, sample.fraction = 0.7, num.threads = as.numeric(cores))
            model_importance_tmp <- as.data.frame(model_tmp$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
            suppressWarnings(model_importance <- merge(model_importance, model_importance_tmp, by = "taxa"))
            
          }
          colnames(model_importance)[2:(nperm + 2)] <- paste0("permutation_", seq(1,nperm + 1))
          model_importance$average <- rowMeans(model_importance[, 2:(nperm + 2)])
          
          ## this else statement does the sample as the above few lines, just for a continous
          ## feature_of_interest...with RF Regression.
        } else {
          model <- ranger::ranger(as.numeric(feature_of_interest) ~ . , data = hData_parent_merge_cor_subset, importance = "impurity_corrected", seed = 42, sample.fraction = 0.7, num.threads = as.numeric(cores))
          model_importance <- as.data.frame(model$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
          for (seed in sample(1:1000, nperm)) {
            model_tmp <- ranger::ranger(as.numeric(feature_of_interest) ~ . , data = hData_parent_merge_cor_subset, importance = "impurity_corrected", seed = seed, sample.fraction = 0.7, num.threads = as.numeric(cores))
            model_importance_tmp <- as.data.frame(model_tmp$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
            suppressWarnings(model_importance <- merge(model_importance, model_importance_tmp, by = "taxa"))
            
          }
          colnames(model_importance)[2:(nperm + 2)] <- paste0("permutation_", seq(1,nperm + 1))
          model_importance$average <- rowMeans(model_importance[, 2:(nperm + 2)])
          
          
        }
        
        ### RF WINNER - PARENT ####
        ## Ok we have the permuted variable.importance. Do, via the average, any
        ## children (sub-species) beat the parents with regards to discriminatory information
        ## with regards to the feature_of_interest? 
        ## If the largest average variable importance is PARENT:
        if ((model_importance %>% dplyr::filter(., taxa %!in% colnames(covariates)) %>% dplyr::arrange(., desc(average)) %>% slice_head(., n = 1) %>% dplyr::pull(., taxa)) == "PARENT") {
          ## drop children from master taxa file taxa_only_split
          mutate_child = children[count]
          taxa_only_split <- taxa_only_split %>% dplyr::mutate(., !!mutate_child := ifelse((.data[[parents[count]]] == parent_trial & !is.na(.data[[parents[count]]]) & na_count == (count - 1)), "drop_dis", .data[[children[count]]]))
          taxa_only_split <- taxa_only_split %>% dplyr::filter(., !grepl(pattern = "drop_dis", x = .data[[children[count]]]))
          assign("taxa_only_split", taxa_only_split, envir = .GlobalEnv)
          
        } else { 
          ### RF WINNER - CHILD ####
          mutate_child = children[count]
          parent_importance <- model_importance$average[model_importance$taxa == "PARENT"]
          
          children_toss <- model_importance %>% dplyr::filter(., taxa %!in% colnames(covariates)) %>% dplyr::filter(., average < parent_importance) %>% dplyr::pull(., taxa)
          children_toss_zero <- model_importance%>% dplyr::filter(., average < 0) %>% dplyr::pull(., taxa) 
          children_toss <- unique(c(children_toss, cor_drop, children_toss_zero))
          
          ## drop parent
          mutate_at_list <- rev(parents)[1:(length(parents) + 1 - count)]
          taxa_only_split <- taxa_only_split %>% 
            dplyr::mutate(., !!mutate_child := ifelse((.data[[parents[count]]] == parent_trial & !is.na(.data[[parents[count]]]) & is.na(.data[[children[count]]])), "drop_dis", .data[[children[count]]]))
          taxa_only_split <- taxa_only_split %>% 
            mutate_at(mutate_at_list, ~ 
                        replace(., .data[[parents[count]]] == parent_trial, NA)) %>%
            dplyr::filter(., !grepl(pattern = "drop_dis", x = .data[[children[count]]]))
          
          ## drop children that didnt win against parent
          taxa_only_split <- taxa_only_split %>% 
            dplyr::mutate(., !!mutate_child := ifelse((.data[[children[count]]] %in% children_toss & is.na(.data[[parents[count]]]) & !is.na(.data[[children[count]]])), "drop_dis", .data[[children[count]]]))
          taxa_only_split <- taxa_only_split %>% 
            dplyr::filter(., !grepl(pattern = "drop_dis", x = .data[[children[count]]]))
          assign("taxa_only_split", taxa_only_split, envir = .GlobalEnv)
          
          
        }
      }
      ### PROGRESS ####
      if (sub_count == length(children_counts[children_counts$count > 1, ][,1] %>% dplyr::pull())) message(paste0("  ", "Done with ", parents[count], " versus ", children[count], ". Moving on to next levels..."))
      sub_count = sub_count + 1
    }
    
    hData <- hData %>%
      dplyr::filter(., clade_name %in% taxa_only_split$clade_name)
    assign("hData", hData, envir = .GlobalEnv)
    count = count + 1
    
  }
  
  ## write no SF files
  hData_output <- hData %>%
    tibble::column_to_rownames(., var = "clade_name") %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(., var = "subject_id") %>%
    janitor::clean_names()
  output_nsf <- merge(metadata, hData_output, by.x = "subject_id", by.y = "subject_id")
  assign("output_nsf", output_nsf, envir = .GlobalEnv)
  assign("output_nsf_count", (NCOL(output_nsf) - 2), envir = .GlobalEnv)
  
  readr::write_delim(file = paste0(tools::file_path_sans_ext(output), "no_sf.txt"), x = output_nsf, delim = "\t")
  readr::write_delim(file = paste0(tools::file_path_sans_ext(output), "_taxa_list_no_sf.txt"), x = taxa_only_split, delim = "\t")
  
  ## assign cor vs rf battle tally to global env
  assign("cor_vs_rf", cor_vs_rf, envir = .GlobalEnv)
 
   ## print results of cor vs rf battle tally
  cat(paste0("\nPercent of parent-child competitions decided through correlation: ", (sum(cor_vs_rf)/ length(cor_vs_rf) * 100), "\n"))
  cat(paste0("Percent of parent-child competitions decided through RF models: ", ((length(cor_vs_rf) - sum(cor_vs_rf))/ length(cor_vs_rf) * 100), "\n"))
}

## super filter covariates =====================================================
## essentially a copy of main super filter function just with added
## merge and consideration of additional features

super_filter_covariates <- function(input = hData, covariates, feature_type = "factor", cores = 4, subsample = opt$subsample, output) {
  
  cat("\n\n", "##################################\n", "Starting Super-filter\n", "##################################\n\n")
  
  ## helper function
  `%!in%` <- Negate(`%in%`)
  
  hData_sf <- input %>%
    tibble::column_to_rownames(., var = "clade_name") %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(., var = "subject_id") %>%
    janitor::clean_names()
  
  hData_sf <- merge(metadata, hData_sf, by.x = "subject_id", by.y = "subject_id")
  hData_sf <- merge(hData_sf, covariates, by = "subject_id")
  covariate_list <- covariates %>% dplyr::select(., -subject_id) %>% colnames()
  output_sf <- hData_sf %>% dplyr::select(., -subject_id) 
  nperm = nperm + 190
  
  pb <- progress_bar$new( format = " Super-filter [:bar] :percent in :elapsed", total = nperm, clear = FALSE, width= 60)
  
  if (feature_type == "factor") {  
    model <- ranger::ranger(as.factor(feature_of_interest) ~ . , data = output_sf, importance = "impurity_corrected", seed = 42, sample.fraction = class_frequencies, replace = TRUE, num.threads = as.numeric(cores))
    model_importance <- as.data.frame(model$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
    for (seed in sample(1:1000, nperm)) {
      
      model_tmp <- ranger::ranger(as.factor(feature_of_interest) ~ . , data = output_sf, importance = "impurity_corrected", seed = seed, sample.fraction = class_frequencies, replace = TRUE, num.threads = as.numeric(cores))
      model_importance_tmp <- as.data.frame(model_tmp$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
      suppressWarnings(model_importance <- merge(model_importance, model_importance_tmp, by = "taxa"))
      
      pb$tick()
    }
    colnames(model_importance)[2:(nperm + 2)] <- paste0("permutation_", seq(1,nperm + 1))
    model_importance$average <- rowMeans(model_importance[, 2:(nperm + 2)])
    model_importance <- model_importance %>% dplyr::relocate(., average)
    assign("model_importance", model_importance, envir = .GlobalEnv)
    
  } else {
    model <- ranger::ranger(as.numeric(feature_of_interest) ~ . , data = output_sf, importance = "impurity_corrected", seed = 42, sample.fraction = as.numeric(subsample), num.threads = as.numeric(cores))
    model_importance <- as.data.frame(model$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
    for (seed in sample(1:1000, nperm)) {
      
      model_tmp <- ranger::ranger(as.numeric(feature_of_interest) ~ . , data = output_sf, importance = "impurity_corrected", seed = seed, sample.fraction = as.numeric(subsample), num.threads = as.numeric(cores))
      model_importance_tmp <- as.data.frame(model_tmp$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
      suppressWarnings(model_importance <- merge(model_importance, model_importance_tmp, by = "taxa"))
      
      pb$tick()
    }
    colnames(model_importance)[2:(nperm + 2)] <- paste0("permutation_", seq(1,nperm + 1))
    model_importance$average <- rowMeans(model_importance[, 2:(nperm + 2)])
    model_importance <- model_importance %>% dplyr::relocate(., average)
    assign("model_importance", model_importance, envir = .GlobalEnv)
  }
  
  
  model_importance_list <- model_importance %>% dplyr::filter(., taxa %!in% covariate_list) %>% dplyr::filter(., average > mean(average)) %>% dplyr::filter(., average > 0) %>% dplyr::pull(., taxa)
  output_sf <- hData_sf %>% dplyr::select(., subject_id, feature_of_interest, all_of(model_importance_list))
  taxa_only_split$clade_name <- janitor::make_clean_names(taxa_only_split$clade_name)
  taxa_only_split <- taxa_only_split %>% dplyr::filter(., clade_name %in% model_importance_list)
  
  assign("output_sf", output_sf, envir = .GlobalEnv)
  assign("taxa_only_split_sf", taxa_only_split, envir = .GlobalEnv)
  
  readr::write_delim(file = output, x = output_sf, delim = "\t")
  readr::write_delim(file = paste0(tools::file_path_sans_ext(output), "_taxa_list.txt"), x = taxa_only_split, delim = "\t")
  
  cat("\n\n Features (no super filter): ", output_nsf_count)
  cat("\n Features (super filter): ", (NCOL(output_sf) - 2), "\n\n")
  
}