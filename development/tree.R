library(data.tree)

hData <- read.delim("/home/docker/example_inputs/microbiome_data.txt", check.names = F)
metadata <- read.delim("/home/docker/example_inputs/metadata.txt", sep = "\t", check.names = FALSE)
hTree <- Node$new("taxaTree", id = 0)

## set random seed if needed
set.seed(42)
nperm <- 10
## helper functions ============================================================

## Negate function ("not in"):
`%!in%` <- Negate(`%in%`)

## suppress warnings
options(warn = -1)


filterPrevalence <- 0.01
filterMeanAbundance <- 0.0001
corrThreshold <- 0.95
trim <- 0.02
feature_of_interest <- "Category"
feature_type <- "factor"
ncores <- 4
sample_fraction <- 1
subject_identifier <- "Sample"

options(width = 150)

metadata <- metadata %>%
  dplyr::select(., subject_identifier, feature_of_interest) %>%
  dplyr::rename(., "feature_of_interest" = feature_of_interest)

# get_descendant_winners goes through the descendants of node, returning a list of all found winners
# maxDepth defines how deep the winner function will go to find a winner
get_descendant_winners <- function(node, max_depth) {
  winners <- list()

  # if maxDepth is zero, this is the bottom
  # return an empty list as no further generations will be considered
  if (max_depth == 0) {
    return(winners)
  }

  for (child in node$children) {
    # if the child is a winner, add the child to list and move to the next child
    if (child$winner) {
      winners <- append(winners, child)
      next
    }

    # otherwise, check the child's children for winners
    winners <- append(winners, get_descendant_winners(child, max_depth - 1))
  }

  return(winners)
}

# build_tree takes a dataframe as input
# rows are the clade names, with levels separated by "|"
# it is expected that there will be a row for every level with summed abundances
# columns are the subjects
# and returns a data.tree structure
build_tree <- function(df, filter_prevalence, filter_mean_abundance) {
  # root node for the tree
  taxa_tree <- data.tree::Node$new("taxaTree", id = 0)

  for (row in seq_len(nrow(df))) {
    # generate a vector of clade levels
    levels <- unlist(strsplit(hData[row, "clade_name"], "\\|"))

    # start at the root node
    node <- taxa_tree

    # iterate the vector of clade level names
    # if the level hasn't been added yet, add it as a new child to the current node
    # otherwise get a reference to the existing node for further iteration
    for (level in levels) {
      potential_node <- node[[level]]
      if (is.null(potential_node)) {
        node <- node$AddChild(level)
      } else {
        node <- potential_node
      }
    }

    # after iterating the levels, node is assigned to the leaf of this row
    # add in the row data and other supporting information

    # a single row dataframe of the abundance data for this row of the df
    node$abundance <- df[row, 2:ncol(df)]
    # indicates if the prevalence filter was passed
    node$passed_prevalence_filter <-
      rowSums(node$abundance != 0) > (NCOL(node$abundance) * filter_prevalence)
    # indicates if the mean abundance filter was passed
    node$passed_mean_abundance_filter <-
      mean(unlist(node$abundance), trim = trim) > filter_mean_abundance
    # defaults to be modified later
    node$winner <- FALSE
    node$highlyCorrelated <- FALSE
    node$lost_rf <- FALSE
    # a unique id corresponding to the row number in the input data
    node$id <- row
  }

  return(taxa_tree)
}

# competes the input node
# this will modify the node in place
# for this node, evaluates correlation and rf against all descendants that have won previous rounds
# modifies the node indicating if it is a winner against those descendants
# OR which of 1:n descendants have out-competed it
compete_node <- function(node, max_depth, corr_threshold, metadata, sample_fraction, ncores) {
  # handle no children, this node is the winner
  if (length(node$children) == 0) {
    node$winner <- TRUE
    return()
  }

  # build dataframe of parent and descendant winners
  # parent is always row 1
  df <- as.data.frame(node$abundance, check.names = FALSE)

  descendant_winners <- get_descendant_winners(node, max_depth)
  # if no descendant winners, th parent is the winner
  # TODO: is this possible? should it be indicated somehow to the end user?
  if (length(descendant_winners) == 0) {
    node$winner <- TRUE
    return()
  }

  # add the descendant's abundance dataframe row to df
  for (descendant in descendant_winners) {
    df <- rbind(df, descendant$abundance)
  }

  # transpose the dataframe to fit the input format for the correlation and ml
  transposed <- as.data.frame(t(df))

  # determine the child ids that are strongly correlated
  correlated_ids <- calculate_correlation(df = transposed, corr_threshold)

  # mark correlated in tree
  # highly correlated descendants are not winners
  for (descendant in descendant_winners) {
    if (descendant$id %in% correlated_ids) {
      descendant$winner <- FALSE
      descendant$highly_correlated <- TRUE
    }
  }

  # if all descendants are correlated, the parent wins
  if (length(descendant_winners) == length(correlated_ids)) {
    node$winner <- TRUE
    return()
  }

  # drop from transposed data all correlated children
  transposed <- transposed %>%
    dplyr::select(., -dplyr::any_of(correlated_ids))

  # run the random forest on the remaining parent + descendants
  # TODO: do the baked-in values need to be modifiable?
  rf_winners <- rf_competition_gpt(df = transposed, metadata, "feature_of_interest", "Sample", "factor", sample_fraction, ncores)

  # mark winners/losers of parent and descendants
  # also mark the non-winners as rf losers
  for (descendant in descendant_winners) {
    if (descendant$id %in% rf_winners) {
      descendant$winner <- TRUE
    } else {
      descendant$winner <- FALSE
      descendant$lost_rf <- TRUE
    }
  }

  return()
}

# competes an entries tree
# takes a data.tree root node as input
# modify_tree: determines if the input tree will be modified in place
# max_depth: determines how deep the descendant competitions will be held
#   defaults to a massive number to allow every descendant
# corr_threshold: the threshold to mark a descendant as highly correlated
# metadata: the metadata associated with the input df that generated the tree
# sample_fraction: TODO (don't know what this is)
# ncores: the number of cores to use when running the random forest
compete_tree <- function(tree, modify_tree = FALSE, max_depth = 1000, corr_threshold, metadata, sample_fraction, ncores) {
  # if not modifying the input tree, create a copy of the tree to perform the competition
  if (!modify_tree) tree <- data.tree::Clone(tree)

  # perform the competition, modifying the tree (which may or may not be a clone of the input)
  tree$Do(
    compete_node,
    max_depth = max_depth,
    corr_threshold = corr_threshold,
    metadata = metadata,
    sample_fraction = sample_fraction,
    ncores = ncores,
    traversal = "post_order"
  )

  # return the tree
  return(tree)
}

calculate_correlation <- function(df, corrThreshold) {
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
      dplyr::mutate(prop = prop.table(n)) %>%
      dplyr::pull(prop)
    ## make the class frequencies a fraction of the entire data to help
    ## prevent overfitting. Janky, but i think this is better than nothing
    ## edit...im gonna leave this in here, but since ranger uses OOB to calc
    ## scores maybe its best to use all the data. Not ideal, but because
    ## sample sizes are usually a little small i think its better to use all
    ## data??
    class_frequencies <- class_frequencies * as.numeric(sample_fraction)
    assign("class_frequencies", class_frequencies, envir = .GlobalEnv)
  } else {
    class_frequencies <- 1
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
    model <- ranger::ranger(as.factor(feature_of_interest) ~ ., data = merged_data, importance = "impurity_corrected", seed = 42, sample.fraction = class_frequencies, replace = TRUE, num.threads = as.numeric(ncores))
    model_importance <- as.data.frame(model$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")

    ## welp, the initial model miiight be correct, but lets permute that process
    ## by looping over some random seeds (10, as set in a very early variable nperm at the top)
    ## and keeping track of the variable.importance. We can then average that and have
    ## a more sure guess whether a Parent or Child is more important with regards to
    ## the feature_of_interest
    for (seed in sample(1:1000, nperm)) {
      model_tmp <- ranger::ranger(as.factor(feature_of_interest) ~ ., data = merged_data, importance = "impurity_corrected", seed = seed, sample.fraction = class_frequencies, replace = TRUE, num.threads = as.numeric(ncores))
      model_importance_tmp <- as.data.frame(model_tmp$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
      suppressWarnings(model_importance <- merge(model_importance, model_importance_tmp, by = "taxa"))
    }
    colnames(model_importance)[2:(nperm + 2)] <- paste0("permutation_", seq(1, nperm + 1))
    model_importance$average <- rowMeans(model_importance[, 2:(nperm + 2)])



    ## this else statement does the sample as the above few lines, just for a continous
    ## feature_of_interest...with RF Regression.
  } else {
    model <- ranger::ranger(as.numeric(feature_of_interest) ~ ., data = merged_data, importance = "impurity_corrected", seed = 42, sample.fraction = sample_fraction, num.threads = as.numeric(ncores))
    model_importance <- as.data.frame(model$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
    for (seed in sample(1:1000, nperm)) {
      model_tmp <- ranger::ranger(as.numeric(feature_of_interest) ~ ., data = merged_data, importance = "impurity_corrected", seed = seed, sample.fraction = sample_fraction, num.threads = as.numeric(ncores))
      model_importance_tmp <- as.data.frame(model_tmp$variable.importance) %>% tibble::rownames_to_column(., var = "taxa")
      suppressWarnings(model_importance <- merge(model_importance, model_importance_tmp, by = "taxa"))
    }
    colnames(model_importance)[2:(nperm + 2)] <- paste0("permutation_", seq(1, nperm + 1))
    model_importance$average <- rowMeans(model_importance[, 2:(nperm + 2)])
  }

  parentColumn <- janitor::make_clean_names(colnames(df)[1])
  if ((model_importance %>% dplyr::arrange(., desc(average)) %>% slice_head(., n = 1) %>% dplyr::pull(., taxa)) == parentColumn) {
    rf_winners <- gsub(pattern = "x", replacement = "", x = parentColumn)
    return(rf_winners)
  } else {
    parent_importance <- model_importance$average[model_importance$taxa == parentColumn]
    children_toss <- model_importance %>%
      dplyr::filter(., average < parent_importance) %>%
      dplyr::pull(., taxa)
    children_toss_zero <- model_importance %>%
      dplyr::filter(., average < 0) %>%
      dplyr::pull(., taxa)
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

hTree <- build_tree(hData, filter_prevalence, filter_mean_abundance)

competed_tree <- compete_tree(
  hTree[["k__Bacteria"]],
  corr_threshold = corr_threshold,
  metadata = metadata,
  sample_fraction = sample_fraction,
  ncores = ncores
)

print(hTree, "winner")


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
