library(data.tree)

hData <- read.delim("/home/docker/example_inputs/microbiome_data.txt", check.names = F)
metadata <- read.delim("/home/docker/example_inputs/metadata.txt", sep = "\t", check.names = FALSE)
hTree <- Node$new("taxaTree", id= 0)

## set random seed if needed
set.seed(42)
nperm = 10
## helper functions ============================================================

## Negate function ("not in"):
`%!in%` <- Negate(`%in%`)

## suppress warnings
options(warn=-1)


filterPrevalence <- 0.01
filterMeanAbundance <- 0.0001
corrThreshold <- 0.95
trim <- 0.02
feature_of_interest <- "Category"
feature_type <- "factor"
ncores <- 4
sample_fraction <- 1
subject_identifier <- "Sample"

options(width=150)

metadata <- metadata %>% dplyr::select(., subject_identifier, feature_of_interest) %>%
  dplyr::rename(., "feature_of_interest" = feature_of_interest)

# go through the descendants of node, returning a list of all found winners
# maxDepth defines how deep the winner function will go to find a winner
getDescendantWinners <- function(node, maxDepth) {
  # if no children or maxDepth is zero, return empty list
  if (length(node$children) == 0 | maxDepth == 0) {
    return(list())
  }
  
  winnerList <- list()
  for (i in 1:length(node$children)) {
   # print(node$name[[i]])
    if (node$children[[i]]$winner) {
      winnerList <- append(winnerList, node$children[[i]])
      next
    }
    winnerList <- append(winnerList, getDescendantWinners(node$children[[i]], maxDepth-1))
  }
  
  return(winnerList)
}

calculateCorrelation <- function(df, corrThreshold) {
  parentColumn <- colnames(df)[1]
  return(
    suppressMessages(corrr::correlate(df)) %>% 
      corrr::focus(., parentColumn) %>% 
      dplyr::filter(., .[[2]] >= as.numeric(corrThreshold)) %>% 
      dplyr::pull(., term)
  )
}

calculateCorrelation_debug <- function(df, corrThreshold) {
  parentColumn <- colnames(df)[1]
  return(
    suppressMessages(corrr::correlate(df)) %>% 
      corrr::focus(., parentColumn)
  )
}

calc_class_frequencies <- function(input, type, feature, sample_fraction) {
  
  if (type == "factor") {
    ## create a vector to deal with class imbalance
    ## this is infomred by: https://github.com/imbs-hl/ranger/issues/167
    ## note in the ranger model, they use sample.fraction and replace = T
    class_frequencies <- input %>% 
      dplyr::count(.data[[feature]]) %>% 
      dplyr::mutate(prop = prop.table(n)) %>% dplyr::pull(prop)
    ## make the class frequencies a fraction of the entire data to help
    ## prevent overfitting. Janky, but i think this is better than nothing
    ## edit...im gonna leave this in here, but since ranger uses OOB to calc
    ## scores maybe its best to use all the data. Not ideal, but because
    ## sample sizes are usually a little small i think its better to use all
    ## data??
    class_frequencies <- class_frequencies * as.numeric(sample_fraction)
    assign("class_frequencies", class_frequencies, envir = .GlobalEnv)
  } else {
    class_frequencies = 1
    assign("class_frequencies", class_frequencies, envir = .GlobalEnv)
  }
  
}

calc_class_frequencies(input = metadata, type = feature_type, feature = "feature_of_interest", sample_fraction = sample_fraction)

rf_competition <- function(df, metadata, feature_of_interest, subject_identifier, feature_type, sample_fraction, ncores) {
  
  merged_data <- merge(df, metadata, by.x = "row.names", by.y = subject_identifier)
  merged_data <- tibble::column_to_rownames(merged_data, var = "Row.names")
  data_colnames <- colnames(merged_data)
  merged_data <- merged_data %>% janitor::clean_names()
  
  if (feature_type == "factor") {
    
    ## create an intial model and keep track of the variable.importance
    ## this will help us decide if the PARENT (species) or the CHILD (sub-speceies)
    ## brings more information to the table with regards to feature_of_interest
    model <- ranger::ranger(as.factor(feature_of_interest) ~ . , data = merged_data, importance = "impurity_corrected", seed = 42, sample.fraction = class_frequencies, replace = TRUE, num.threads = as.numeric(ncores))
    model_importance <- as.data.frame(model$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
    
    ## welp, the initial model miiight be correct, but lets permute that process
    ## by looping over some random seeds (10, as set in a very early variable nperm at the top)
    ## and keeping track of the variable.importance. We can then average that and have
    ## a more sure guess whether a Parent or Child is more important with regards to 
    ## the feature_of_interest
    for (seed in sample(1:1000, nperm)) {
      model_tmp <- ranger::ranger(as.factor(feature_of_interest) ~ . , data = merged_data, importance = "impurity_corrected", seed = seed, sample.fraction = class_frequencies, replace = TRUE, num.threads = as.numeric(ncores))
      model_importance_tmp <- as.data.frame(model_tmp$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
      suppressWarnings(model_importance <- merge(model_importance, model_importance_tmp, by = "taxa"))
      
    }
    colnames(model_importance)[2:(nperm + 2)] <- paste0("permutation_", seq(1,nperm + 1))
    model_importance$average <- rowMeans(model_importance[, 2:(nperm + 2)])
    
    
    
    ## this else statement does the sample as the above few lines, just for a continous
    ## feature_of_interest...with RF Regression.
  } else {
    model <- ranger::ranger(as.numeric(feature_of_interest) ~ . , data = merged_data, importance = "impurity_corrected", seed = 42, sample.fraction = sample_fraction, num.threads = as.numeric(ncores))
    model_importance <- as.data.frame(model$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
    for (seed in sample(1:1000, nperm)) {
      model_tmp <- ranger::ranger(as.numeric(feature_of_interest) ~ . , data = merged_data, importance = "impurity_corrected", seed = seed, sample.fraction = sample_fraction, num.threads = as.numeric(ncores))
      model_importance_tmp <- as.data.frame(model_tmp$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
      suppressWarnings(model_importance <- merge(model_importance, model_importance_tmp, by = "taxa"))
      
    }
    colnames(model_importance)[2:(nperm + 2)] <- paste0("permutation_", seq(1,nperm + 1))
    model_importance$average <- rowMeans(model_importance[, 2:(nperm + 2)])
    
  }
  
  parentColumn <- janitor::make_clean_names(colnames(df)[1])
  if ((model_importance %>% dplyr::arrange(., desc(average)) %>% slice_head(., n = 1) %>% dplyr::pull(., taxa)) == parentColumn) { 
    rf_winners <- gsub(pattern = "x", replacement = "", x = parentColumn)
    return(rf_winners) 
    
  } else {
    parent_importance <- model_importance$average[model_importance$taxa == parentColumn]
    children_toss <- model_importance %>% dplyr::filter(., average < parent_importance) %>% dplyr::pull(., taxa)
    children_toss_zero <- model_importance %>% dplyr::filter(., average < 0) %>% dplyr::pull(., taxa) 
    children_toss <- unique(c(children_toss, children_toss_zero))
    children_winners <- model_importance$taxa[model_importance$taxa %!in% c(children_toss, parentColumn)]
    rf_winners <- gsub(pattern = "x", replacement = "", x = children_winners)
    return(rf_winners)
  }
}

rf_competition_gpt <- function(df, metadata, feature_of_interest = "feature_of_interest", subject_identifier, feature_type, sample_fraction = class_frequencies, ncores, nperm = 10) {
  
  merged_data <- merge(df, metadata, by.x = "row.names", by.y = subject_identifier)
  merged_data <- tibble::column_to_rownames(merged_data, var = "Row.names")
  data_colnames <- colnames(merged_data)
  merged_data <- merged_data %>% janitor::clean_names()
  
  if (feature_type == "factor") {
    response_formula <- as.formula(paste("as.factor(", feature_of_interest, ") ~ .", sep = ""))
  } else {
    response_formula <- as.formula(paste("as.numeric(", feature_of_interest, ") ~ .", sep = ""))
  }
  
  run_ranger <- function(seed) {
    ranger::ranger(response_formula, data = merged_data, importance = "impurity_corrected", seed = seed, sample.fraction = sample_fraction, replace = TRUE, num.threads = as.numeric(ncores))$variable.importance %>%
      as.data.frame() %>%
      dplyr::rename(., "importance" = ".") %>%
      tibble::rownames_to_column(var = "taxa")
  }
  
  model_importance <- purrr::map_df(sample(1:1000, nperm), run_ranger) %>%
    dplyr::group_by(taxa) %>%
    dplyr::summarise(., average = mean(importance))

  parentColumn <- janitor::make_clean_names(colnames(df)[1])
  
  if ((model_importance %>% arrange(desc(average)) %>% pull(taxa))[1] == parentColumn) { 
    return(gsub(pattern = "x", replacement = "", x = parentColumn))
  } else {
    parent_importance <- model_importance$average[model_importance$taxa == parentColumn]
    children_toss <- model_importance %>%
      dplyr::filter(average < parent_importance | average < 0) %>%
      dplyr::pull(taxa)
    children_winners <- model_importance %>%
      dplyr::filter(!taxa %in% c(children_toss, parentColumn)) %>%
      dplyr::pull(taxa)
    return(gsub(pattern = "x", replacement = "", x = children_winners))
  }
}

for (row in 1:nrow(hData)) {
  levels <- unlist(strsplit(hData[row, "clade_name"], "\\|"))
  
  node <- hTree
  for (i in 1:length(levels)) {
    exists <- "FALSE"
    if (length(node$children) < 1) {
      node <- node$AddChild(levels[i])
      next
    }
    
    potentialNode <- node[[levels[i]]]
    if (is.null(potentialNode$name) == "TRUE") {
      node <- node$AddChild(levels[i])
    } else {
      node <- potentialNode
    }
  }
  
  node$abundance <- hData[row,2:ncol(hData)]
  node$passedPrevalenceFilter <- rowSums(node$abundance != 0) > (NCOL(node$abundance)*filterPrevalence)
  node$passedMeanAbundanceFilter <- mean(unlist(node$abundance), trim = trim) > filterMeanAbundance
  node$winner <- FALSE    
  node$highlyCorrelated = FALSE
  node$rf_winners = FALSE
  node$id <- row
}

hTree[["k__Bacteria"]]$Do(function(node) node$competed <- {
  # handle no children
  if (length(node$children) == 0) {
    node$winner <- TRUE
    return(TRUE)
  }
  
  # build dataframe of parent and children
  # parent is row 1
  df <- as.data.frame(node$abundance, check.names = FALSE)
  
  # handle descendant winners
  # take in depth limit
  descendantWinners <- getDescendantWinners(node, 6)
  if (length(descendantWinners) == 0) {
    node$winner <- TRUE
    return(TRUE)
  }
  for (i in 1:length(descendantWinners)) {
    df <- rbind(df, descendantWinners[[i]]$abundance)
  }
  
  # transpose
  transposed <- as.data.frame(t(df))
  
  # standalone func
  # inputs: transposed dataframe, correlation threshold
  # mark correlated values
  # outputs: vector of correlated ids
  correlatedIDs <- calculateCorrelation(df = transposed, corrThreshold)
  correlated_df_debug <- calculateCorrelation_debug(df = transposed, corrThreshold)
  
  # mark correlated in tree
  for (i in 1:length(descendantWinners)) {
    if (descendantWinners[[i]]$id %in% correlatedIDs) {
      descendantWinners[[i]]$winner <- FALSE
      descendantWinners[[i]]$highlyCorrelated <- TRUE
    }
  }
  # handle all children correlated, parent winner
  if (length(descendantWinners) == length(correlatedIDs)) {
    node$winner <- TRUE
    return(TRUE)
  }
  # drop from transposed data all correlated children
  transposed <- transposed %>%
    dplyr::select(., -dplyr::any_of(correlatedIDs))
  # merge metadata??
  
  # standalone func
  # inputs: transposed dataframe, metadata, ...
  # random forest step
  # select factor
  # run model
  # outputs: vector of winner ids
  if (length(colnames(transposed)) > 1) {
    rf_winners <- rf_competition_gpt(df = transposed, metadata, "feature_of_interest", "Sample", "factor", sample_fraction, ncores)
  }
  # 
  # mark winners/losers of parent and descendants
  for (i in 1:length(descendantWinners)) {
    if (descendantWinners[[i]]$id %in% rf_winners) {
      descendantWinners[[i]]$winner <- FALSE
      descendantWinners[[i]]$rf_winners <- TRUE
    }
  }
  
}, traversal = "post-order")

flatten_tree_with_metadata <- function(node) {
  df <- data.frame(
    name = node$name,
    depth = node$level,
    pathString = node$pathString,
    rf_win = node$rf_winners,
    highly_cor = node$highlyCorrelated,
    passed_prevelance = node$passedPrevalenceFilter,
    passed_abundance = node$passedMeanAbundanceFilter,
    abundance = node$abundance,
    stringsAsFactors = FALSE
  )
  
  if (length(node$children) > 0) {
    children_df <- do.call(rbind, lapply(node$children, flatten_tree_with_metadata))
    df <- rbind(df, children_df)
  }
  
  return(df)
}

# Flatten the tree and tree decisions
flattened_df <- flatten_tree_with_metadata(hTree[["k__Bacteria"]])
flattened_df <- flattened_df %>% filter(., rf_win == TRUE, passed_prevelance == TRUE, passed_abundance == TRUE, highly_cor == FALSE)





hello <- as.Node(flattened_df)
View(as.data.frame(hello))

