## HFE FUNCTIONS
## v2

library(janitor, quietly = T, verbose = F, warn.conflicts = F)
library(tidyr, quietly = T, verbose = F, warn.conflicts = F)
library(tibble, quietly = T, verbose = F, warn.conflicts = F)
library(ggplot2, quietly = T, verbose = F, warn.conflicts = F)
library(ggsci, quietly = T, verbose = F, warn.conflicts = F)
library(progress, quietly = T, verbose = F, warn.conflicts = F)
library(data.tree, quietly = T, verbose = F, warn.conflicts = F)
library(dplyr, quietly = T, verbose = F, warn.conflicts = F)
library(corrr, quietly = T, verbose = F, warn.conflicts = F)
library(tibble, quietly = T, verbose = F, warn.conflicts = F)
library(purrr, quietly = T, verbose = F, warn.conflicts = F)
library(ranger, quietly = T, verbose = F, warn.conflicts = F)
library(vroom, quietly = T, verbose = F, warn.conflicts = F)
library(tidyselect, quietly = T, verbose = F, warn.conflicts = F)
library(data.table, quietly = T, verbose = F, warn.conflicts = F)
library(dtplyr, quietly = T, verbose = F, warn.conflicts = F)
library(lineprof, quietly = T, verbose = F, warn.conflicts = F)


## set random seed if needed
set.seed(Sys.time())
nperm <- 20 # permute the random forest this many times
trim <- 0.02 # trim outliers from mean feature abundance calc

## helper functions ============================================================

## Negate function ("not in"):
`%!in%` <- Negate(`%in%`)

## suppress warnings
options(warn = -1)

## read in metadata  ===========================================================
read_in_metadata <- function(input, subject_identifier, label) {
  if (strsplit(basename(input), split = "\\.")[[1]][2] %in% c("tsv","txt")) {
    delim = "\t"
  } else {
    delim = ","
  }
  metadata <- suppressMessages(readr::read_delim(file = input, delim = delim)) %>%
    dplyr::select(., subject_identifier, label) %>%
    dplyr::rename(., "subject_id" = subject_identifier) %>%
    rename(., "feature_of_interest" = label) %>%
    tidyr::drop_na()
  ## this is an effort to clean names for the RF later, which will complain big time
  ## if there are symbols or just numbers in your subject IDs
  cat("\n\nIF your subject_identifiers are just numerics OR have weird symbols\nin them, there's a chance this program crashes. The RF engine,\nranger, is a little picky and very stubborn with names.\n")
  cat("TaxaHFE will do its best to handle the names initially...\n")
  cat("(we <3 ranger though)\n")
  metadata$subject_id <- metadata$subject_id %>% janitor::make_clean_names(use_make_names = F)
  return(metadata)
}

## read in microbiome data =====================================================
read_in_microbiome <- function(input, meta = metadata, abundance, format_metaphlan, cores = opt$ncores) {
  
  ## read in txt, tsv, or csv microbiome data
  if (strsplit(basename(input), split = "\\.")[[1]][2] %in% c("tsv","txt")) {
    delim = "\t"
  } else {
    delim = ","
  }
  
  hData <- suppressMessages(vroom::vroom(file = input, delim = delim, skip = 0, 
                                         .name_repair = "minimal", 
                                         num_threads = as.numeric(cores)) %>% 
                              dplyr::select(., -any_of(c("NCBI_tax_id", 
                                                         "clade_taxid"))) %>%
                              janitor::clean_names(use_make_names = F))
  ## only select columns that are in metadata file, reduce computation
  hData <- hData %>% 
    dplyr::select(., dplyr::any_of(c("clade_name", metadata$subject_id)))
  
  ## check and make sure clade_name is the first column in hData
  if ("clade_name" %!in% colnames(hData)) {
    stop("column clade_name not found in input data")
  }
  if (colnames(hData)[1] != "clade_name") {
    hData <- hData %>% dplyr::relocate(., "clade_name")
  }
  
  ## Applying initial abundance cutoffs. This will vastly shrink the dataset usually
  cat("\nApplying initial abundance filters...\n")
  ## count number of splits (+1) in hierarchical data by "|" symbol
  
  #3 get the base level sums per sample
  num_levels <- max(stringr::str_count(hData$clade_name, "\\|"))
  levels <- paste0("L", seq(1:(num_levels + 1)))
  
  if (format_metaphlan == TRUE) {
    hData_total_abundance <- hData %>% 
      dplyr::filter(., !grepl("\\|", clade_name)) %>% 
      tibble::column_to_rownames(., var = "clade_name") %>% 
      summarise_all(sum) %>% t() %>% as.data.frame() %>% rename(., "row_sums" = "V1")
  } else {
  hData_total_abundance <- hData %>% 
    tidyr::separate(., col = "clade_name", into = levels,
                    extra = "drop", sep = "\\|") %>%
    dplyr::select(., tidyselect::where(is.numeric)) %>%
    summarise_all(sum) %>% t() %>% as.data.frame() %>% rename(., "row_sums" = "V1")
  }
  
  ## divide those base level sums by all the features in the dataset
  ## this is a sample by sample relative abundance
  hData_rel_abundance <- sweep((hData %>% tibble::column_to_rownames(., var = "clade_name") 
                                %>% as.matrix), 2, hData_total_abundance$row_sums, "/")
  
  ## filter features with mean abundance that is above the threshold
  high_abundant_taxa <- hData_rel_abundance %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(., var = "clade_name") %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(., resistant_row_means = mean(dplyr::c_across(2:NCOL(hData_rel_abundance)), trim = trim, na.rm = T)) %>% 
    dplyr::filter(., resistant_row_means > as.numeric(abundance)) %>% 
    dplyr::pull(., clade_name)
  
  ## select those features from the original dataset, back to original counts
  hData <- hData %>% dplyr::filter(., clade_name %in% high_abundant_taxa)
  
  ## write input to file
  return(as.data.frame(hData))
}

## write separate files to test summary levels =================================
# write summarized abundance files for each level except taxa_tree
write_summary_files <- function(input, metadata, output) {
  
  max_levels <- max(input[["depth"]])
  ## start at 2 to ignore taxa_tree depth (meaningless node)
  levels <- c(1:max_levels)
  
  ## split raw data by pipe symbol into number of expected parts
  
  count = 1
  for (i in seq(levels)) {
    if (i == 1) {
      next
    }
    ## select different levels and write them to file
    file_summary <- input %>%
      dplyr::filter(., depth == i & passed_prevelance == TRUE & passed_abundance == TRUE) %>%
      dplyr::select(., name, 10:dplyr::last_col()) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames(., var = "name") %>%
      t() %>%
      as.data.frame() %>%
      tibble::rownames_to_column(., var = "subject_id") %>%
      janitor::clean_names()
    
    file_merge <- merge(metadata, file_summary, by = "subject_id")
    
    filename <- paste0("_level_", count, ".csv")
    readr::write_delim(x = file_merge, file = paste0(tools::file_path_sans_ext(output), filename), delim = ",")
    count = count + 1
  }
  
}

## write files for old_HFE =====================================================
# write old files for the Oudah program
write_old_hfe <- function(input, output) {
  
  max_levels <- max(input[["depth"]])
  ## start at 2 to ignore taxa_tree depth (meaningless node)
  levels <- c(1:max_levels)
  
  ## split raw data by pathString backslash into number of expected parts

  taxonomy <- input %>% 
    dplyr::filter(., depth == max_levels & passed_prevelance == TRUE & passed_abundance == TRUE) %>%
    dplyr::select(., pathString) %>%
    tidyr::separate(., col = "pathString", into = c(unlist(paste0("L", c(1:max_levels)))),
                    extra = "drop", sep = "\\/") %>%
    dplyr::select(., -L1) %>%
    tibble::remove_rownames()

  abundance <- input %>%
    dplyr::filter(., depth == max_levels & passed_prevelance == TRUE & passed_abundance == TRUE) %>%
    dplyr::select(., name, 10:dplyr::last_col()) %>%
    tibble::remove_rownames()
    
  taxonomy <- taxonomy[taxonomy[,ncol(taxonomy)] %in% abundance$name, ]
  input_taxa_merge <- merge(taxonomy, abundance, by.x = paste0("L", max_levels), by.y = "name")
  
  input_taxa_merge$index <- (1001:(NROW(input_taxa_merge) + 1000))
  input_taxa_merge$L2 <- "k__Bacteria"
  
  input_taxa_merge <- input_taxa_merge %>%
    dplyr::relocate(., index, dplyr::any_of(c(unlist(paste0("L", c(1:max_levels))))))
  
  readr::write_delim(x = input_taxa_merge[1:(max_levels)], file = paste0(tools::file_path_sans_ext(output), "_old_hfe_taxa.txt"), col_names = FALSE, delim = "\t")
  readr::write_delim(x = input_taxa_merge %>% dplyr::select(., 1,(max_levels+1):dplyr::last_col()), file = paste0(tools::file_path_sans_ext(output), "_old_hfe_otu.txt"), col_names = FALSE, delim = "\t")
  
  metadata_order <- colnames(input_taxa_merge[,9:NCOL(input_taxa_merge)])
  
  metadata_list <- metadata %>% dplyr::arrange(match(subject_id, metadata_order)) %>%
    pull(., feature_of_interest)
  
  metadata_list <- as.data.frame(c("label", metadata_list))
  readr::write_delim(x = as.data.frame(t(metadata_list)), file = paste0(tools::file_path_sans_ext(output), "_old_hfe_label.txt"), col_names = FALSE, delim = "\t")
  
}

## get descenant winners =======================================================
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

## set the initial values for a leaf ===========================================
# modifies the passed-in node
# node: current leaf node
# row_num: unique row number from original df
# row_df: single row dataframe of abundances for this leaf node
# filter_prevalence: passed in filter cutoff for prevalence percentage
# filter_mean_abundance: passed in filter cutoff for mean abundance
initial_leaf_values <- function(node, row_num, row_df, filter_prevalence, filter_mean_abundance) {
  # a single row dataframe of the abundance data for this row of the df
  node$abundance <- row_df
  # indicates if the prevalence filter was passed
  node$passed_prevalence_filter <-
    rowSums(node$abundance != 0) > (NCOL(node$abundance) * as.numeric(filter_prevalence))
  # indicates if the mean abundance filter was passed
  node$passed_mean_abundance_filter <-
    mean(unlist(node$abundance), trim = trim) > filter_mean_abundance
  # defaults to be modified later
  node$winner <- FALSE
  node$highly_correlated <- FALSE
  node$lost_rf <- FALSE
  # list of outcomes, one for every time the node participates in compete_node
  node$outcomes <- list()
  # a unique id corresponding to the row number in the input data
  node$id <- row_num
}

## fix unpopulated node ========================================================
# generates abundance, and other values, for clade nodes that were missing a row in the input data
# uses the sum of the child abundances for the abundance data
# if no child abundances are available, uses all zeros
# node: current node, passed in by the tree traversal function
# zeros_df: a single row df of zeros that matches the original dataframe
# next_row_id: the next id after the last row in the original dataframe

fix_unpopulated_node <- function(node, zeros_df, next_row_id, filter_prevalence, filter_mean_abundance) {
  # Ignore nodes with abundance or no children
  if (!is.null(node$abundance)) {
    return()
  }
  print(paste0("Adding nodes to populate tree for: ", node$name))
  
  # Collect non-null child abundances
  child_abundances <- lapply(node$children, function(child) child$abundance)
  child_abundances <- Filter(function(abundance) !is.null(abundance), child_abundances)
  
  if (length(child_abundances) == 0) {
    # No child abundances, use a row of zeros
    new_row <- as.data.table(zeros_df)[, lapply(.SD, sum)]
  } else {
    # Combine child abundances and calculate row sums
    combined_abundances <- rbindlist(child_abundances)
    new_row <- combined_abundances[, lapply(.SD, sum)]
  }
  
  # Assign unique id to row name
  rownames(new_row) <- next_row_id

  # Populate values now that the summed abundance exists
  initial_leaf_values(node, next_row_id, new_row, filter_prevalence, filter_mean_abundance)
}

## build tree from df ==========================================================
# build_tree takes a dataframe as input
# rows are the clade names, with levels separated by "|"
# it is expected that there will be a row for every level with summed abundances
# columns are the subjects
# and returns a data.tree structure
build_tree <- function(df, filter_prevalence, filter_mean_abundance) {
  # root node for the tree
  taxa_tree <- data.tree::Node$new("taxaTree", id = 0)
  
  ## progress bar
  pb <- progress_bar$new( format = " Adding nodes to tree [:bar] :percent in :elapsed", total = nrow(df), clear = FALSE, width = 60)
  
  for (row in seq_len(nrow(df))) {
    ## progress bar
    pb$tick()
    
    # generate a vector of clade levels
    levels <- unlist(strsplit(df[row, "clade_name"], "\\|"))
    
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
    initial_leaf_values(node, row, df[row, 2:ncol(df)], filter_prevalence, filter_mean_abundance)
    
  }
  
  # now that the tree is built, handle unpopulated leaves with the fix_unpopulated_node
  # generate a single row zeros df for a default abundance in the case of no child data to sum
  # matches other abundance by no including clade_name column
  zeros_df <- df[1, 2:ncol(df)]
  zeros_df[zeros_df != 0] <- 0
  
  # start the unique id counter at 1 greater than the original df size
  next_row_id <- nrow(df) + 1
  
  # traverse the tree and fix the unpopulated nodes
  taxa_tree$Do(function(node) {
    fix_unpopulated_node(node, zeros_df, next_row_id, filter_prevalence, filter_mean_abundance)
    # this loop handles that by incrementing next_row_id for every node, ensuring a unique id if needed
    # they are NOT guaranteed to be sequential since the current node may or may not need it
    # <<- ensures that we assign to the next_row_id var outside this closure loop
    next_row_id <<- next_row_id + 1
  }, traversal = "post-order")

  return(taxa_tree)
}

## main compete node function ==================================================
# competes the input node
# this will modify the node in place
# for this node, evaluates correlation and rf against all descendants that have won previous rounds
# modifies the node indicating if it is a winner against those descendants
# OR which of 1:n descendants are winners
compete_node <- function(node, lowest_level, max_depth, corr_threshold, metadata, sample_fraction, ncores, feature_type, nperm) {
  # skip anything lower than the lowest level (exclusive)
  if (node$level < lowest_level) {
    return()
  }
  
  ## do not consider children that do not pass abundance and prevalence filters
  if (!node$passed_prevalence_filter || !node$passed_mean_abundance_filter) {
    node$outcomes <- append(node$outcomes, "loss: did not pass filters")
    return()
  }
  
  # handle no children, this node is the winner
  if (length(node$children) == 0) {
    node$outcomes <- append(node$outcomes, "win: no children")
    node$winner <- TRUE
    return()
  }
  
  # build dataframe of parent and descendant winners
  # parent is always row 1
  df <- as.data.frame(node$abundance, check.names = FALSE)
  
  descendant_winners <- get_descendant_winners(node, max_depth)
  # if no descendant winners, the parent is the winner
  # TODO: is this possible? should it be indicated somehow to the end user?
  if (length(descendant_winners) == 0) {
    node$outcomes <- append(node$outcomes, "win: no descendant winners")
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
  not_correlated_descendant_winners <- list()
  for (descendant in descendant_winners) {
    if (descendant$id %in% correlated_ids) {
      descendant$outcomes <- append(descendant$outcomes, sprintf("loss: highly correlated to ancestor %s", node$name))
      descendant$winner <- FALSE
      descendant$highly_correlated <- TRUE
    } else {
      not_correlated_descendant_winners <- append(not_correlated_descendant_winners, descendant)
    }
  }
  
  # if all descendants are correlated, the parent wins
  if (length(descendant_winners) == length(correlated_ids)) {
    node$outcomes <- append(node$outcomes, sprintf(
      "win: all descendant winners highly correlated, %s",
      paste(
        sapply(
          descendant_winners,
          function(node) node$name
        ),
        collapse = ", "
      )
    ))
    node$winner <- TRUE
    return()
  }
  
  # drop from transposed data all correlated children
  transposed <- transposed %>%
    dplyr::select(., -dplyr::any_of(correlated_ids))
  
  # run the random forest on the remaining parent + descendants
  # TODO: do the baked-in values need to be modifiable?
  rf_winners <- rf_competition(transposed, metadata, "feature_of_interest", "subject_id", feature_type, sample_fraction, ncores, nperm)
  
  # generate winner and loser name lists from the competitors (parent and non-correlated descendants), using the outcome
  # build ahead of time so that a summary can be provided in outcomes
  # this can be sped up and done in a single loop if outcomes are not needed
  competitors <- append(not_correlated_descendant_winners, node)
  winner_names <- list()
  loser_names <- list()
  for (competitor in competitors) {
    if (competitor$id %in% rf_winners) {
      winner_names <- append(winner_names, competitor$name)
    } else {
      loser_names <- append(loser_names, competitor$name)
    }
  }
  outcome_str <- sprintf("winners: %s; losers: %s", paste(winner_names, collapse = ","), paste(loser_names, collapse = ","))
  
  # now actually mark the results, including outcome string
  # mark winners/losers of parent and descendants
  # also mark the non-winners as rf losers
  for (competitor in competitors) {
    if (competitor$id %in% rf_winners) {
      competitor$outcomes <- append(competitor$outcomes, sprintf("win: rf winner, %s", outcome_str))
      competitor$winner <- TRUE
    } else {
      competitor$outcomes <- append(competitor$outcomes, sprintf("loss: rf loser, %s", outcome_str))
      competitor$winner <- FALSE
      competitor$lost_rf <- TRUE
    }
  }
  
  return()
}

## compete all winners (final RF) ==============================================
# compete all winners, updating the tree in the process
# TODO: combine the overlaps in this code with the code in the function above
compete_all_winners <- function(tree, metadata, sample_fraction, ncores) {
  # all vs all competition with winners
  # skipped rows have winner = FALSE so won't appear in this list
  competitors <- get_descendant_winners(tree, tree$height)
  if (length(competitors) == 0) {
    return()
  }
  
  # all winners into a transposed data frame
  df <- data.frame()
  for (winner in competitors) {
    df <- rbind(df, winner$abundance)
  }
  transposed <- as.data.frame(t(df))
  
  # compete here
  # return list of winner ids
  # right now does nothing, all are marked as winners
  rf_winners <- lapply(competitors, function(node) node$id)
  
  # TODO: so much duplication below
  
  # generate winner and loser name lists from the competitors (parent and non-correlated descendants), using the outcome
  # build ahead of time so that a summary can be provided in outcomes
  # this can be sped up and done in a single loop if outcomes are not needed
  winner_names <- list()
  loser_names <- list()
  for (competitor in competitors) {
    if (competitor$id %in% rf_winners) {
      winner_names <- append(winner_names, competitor$name)
    } else {
      loser_names <- append(loser_names, competitor$name)
    }
  }
  outcome_str <- sprintf("winners: %s; losers: %s", paste(winner_names, collapse = ","), paste(loser_names, collapse = ","))
  
  # now actually mark the results, including outcome string
  # mark winners/losers of parent and descendants
  # also mark the non-winners as rf losers
  for (competitor in competitors) {
    if (competitor$id %in% rf_winners) {
      competitor$outcomes <- append(competitor$outcomes, sprintf("win: rf winner, %s", outcome_str))
      competitor$winner <- TRUE
    } else {
      competitor$outcomes <- append(competitor$outcomes, sprintf("loss: rf loser, %s", outcome_str))
      competitor$winner <- FALSE
      competitor$lost_rf <- TRUE
    }
  }
}

## wrapper function to loop through tree nodes =================================
# competes an entries tree
# takes a data.tree root node as input
# modify_tree: determines if the input tree will be modified in place
# lowest_level: lowest level of the tree to consider during operations
#   this level will be compared in a final all-vs-all random forest after the main tree competition
#   defaults to skipping the lowest level (all abundances)
# max_depth: determines how deep the descendant competitions will be held
#   defaults to a massive number to allow every descendant
# corr_threshold: the threshold to mark a descendant as highly correlated
# metadata: the metadata associated with the input df that generated the tree
# sample_fraction: fraction of data to use in rf to help prevent data leakage
# ncores: the number of cores to use when running the random forest
compete_tree <- function(tree, modify_tree = TRUE, lowest_level = 2, max_depth = 1000, corr_threshold, metadata, sample_fraction, ncores, feature_type, nperm) {
  # if not modifying the input tree, create a copy of the tree to perform the competition
  if (!modify_tree) tree <- data.tree::Clone(tree)
  
  # perform the competition, modifying the tree (which may or may not be a clone of the input)
  tree$Do(
    compete_node,
    lowest_level = lowest_level,
    max_depth = max_depth,
    corr_threshold = corr_threshold,
    metadata = metadata,
    sample_fraction = sample_fraction,
    ncores = ncores,
    feature_type = feature_type,
    nperm = nperm,
    traversal = "post-order"
  )
  
  # compete all winners
  compete_all_winners(tree, metadata, sample_fraction, ncores)
  
  # return the tree
  return(tree)
}

## calculate node correlation ==================================================
# calculate correlation between nodes and the remaining children
calculate_correlation <- function(df, corrThreshold) {
  parentColumn <- colnames(df)[1]
  return(
    suppressMessages(corrr::correlate(df)) %>%
      corrr::focus(., parentColumn) %>%
      dplyr::filter(., .[[2]] >= as.numeric(corrThreshold)) %>%
      dplyr::pull(., term)
  )
}

## calculate class frequencies =================================================
# function to calculate frequencies of factor classes from metadata
# and pass that into the ranger function, to help with class imbalance. These
# class frequencies are modified by the subsample flag (a subsample of the
# class frequencies).
# If metadata variable is not a factor, returns the subsample flag
calc_class_frequencies <- function(input, feature_type, feature = "feature_of_interest", sample_fraction) {
  if (feature_type == "factor") {
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
    return(class_frequencies)
  } else {
    return(as.numeric(sample_fraction))
  }
}

## rf competition function =====================================================
rf_competition <- function(df, metadata, feature_of_interest = "feature_of_interest", subject_identifier = "subject_id", feature_type = feature_type, sample_fraction = calc_class_frequencies(), ncores = ncores, nperm = nperm) {
  # merge node abundance + children abundance with metadata
  merged_data <- merge(df, metadata, by.x = "row.names", by.y = "subject_id")
  # clean node names so ranger doesnt throw an error
  merged_data <- tibble::column_to_rownames(merged_data, var = "Row.names")
  data_colnames <- colnames(merged_data)
  merged_data <- merged_data %>% janitor::clean_names()
  
  # determine if rf regression or classification should be run
  if (feature_type == "factor") {
    response_formula <- as.formula(paste("as.factor(", feature_of_interest, ") ~ .", sep = ""))
  } else {
    response_formula <- as.formula(paste("as.numeric(", feature_of_interest, ") ~ .", sep = ""))
  }
  
  # run ranger, setting parameters such as
  # random seed
  # sample.fraction (acquired from class frequencies function)
  # num.threads number of threads to five ranger
  run_ranger <- function(seed) {
    ranger::ranger(response_formula, data = merged_data, importance = "impurity_corrected", seed = seed, sample.fraction = sample_fraction, replace = TRUE, num.threads = as.numeric(ncores))$variable.importance %>%
      as.data.frame() %>%
      dplyr::rename(., "importance" = ".") %>%
      tibble::rownames_to_column(var = "taxa")
  }
  
  # run the above function across 10 random seeds and average the vip scores
  model_importance <- purrr::map_df(sample(1:1000, nperm), run_ranger) %>%
    dplyr::group_by(taxa) %>%
    dplyr::summarise(., average = mean(importance))
  
  # specify the parent column, which is the vip score to beat
  parentColumn <- janitor::make_clean_names(colnames(df)[1])
  
  # if top vip score is the parent, parent wins, else grab the children who 
  # beat the parent's VIP score
  if ((model_importance %>% arrange(desc(average)) %>% pull(taxa))[1] == parentColumn) {
    return(gsub(pattern = "x", replacement = "", x = parentColumn))
  } else {
    parent_importance <- model_importance$average[model_importance$taxa == parentColumn]
    children_toss <- model_importance %>%
      dplyr::filter(average < parent_importance) %>%
      dplyr::pull(taxa)
    children_winners <- model_importance %>%
      dplyr::filter(!taxa %in% c(children_toss, parentColumn)) %>%
      dplyr::pull(taxa)
    return(gsub(pattern = "x", replacement = "", x = children_winners))
  }
}

## Flatten tree to data frame ==================================================
## exports tree as dataframe with tons of info on how the competition went
flatten_tree_with_metadata <- function(node) {
  df <- data.frame(
    name = node$name,
    depth = node$level,
    pathString = node$pathString,
    outcomes = paste(node$outcomes, collapse = "|\n"),
    winner = node$winner,
    rf_loss = node$lost_rf,
    highly_cor = node$highly_correlated,
    passed_prevelance = node$passed_prevalence_filter,
    passed_abundance = node$passed_mean_abundance_filter,
    abundance = node$abundance,
    stringsAsFactors = FALSE
  )

  if (length(node$children) > 0) {
    children_df <- do.call(rbind, lapply(node$children, flatten_tree_with_metadata))
    df <- rbind(df, children_df)
  }
  return(df)
}

## super filter ================================================================

rf_competition_sf <- function(df, metadata, 
                              feature_of_interest = "feature_of_interest", 
                              subject_identifier = "subject_id", 
                              feature_type = feature_type, 
                              sample_fraction = calc_class_frequencies(), 
                              ncores = ncores, nperm = nperm,
                              output) {

  # determine if rf regression or classification should be run
  if (feature_type == "factor") {
    response_formula <- as.formula(paste("as.factor(", feature_of_interest, ") ~ .", sep = ""))
  } else {
    response_formula <- as.formula(paste("as.numeric(", feature_of_interest, ") ~ .", sep = ""))
  }
  
  # run ranger, setting parameters such as
  # random seed
  # sample.fraction (acquired from class frequencies function)
  # num.threads number of threads to five ranger
  run_ranger <- function(seed) {
    ranger::ranger(response_formula, data = df, importance = "impurity_corrected", seed = seed, sample.fraction = sample_fraction, replace = TRUE, num.threads = as.numeric(ncores))$variable.importance %>%
      as.data.frame() %>%
      dplyr::rename(., "importance" = ".") %>%
      tibble::rownames_to_column(var = "taxa")
  }
  
  # run the above function across 10 random seeds and average the vip scores
  model_importance <- purrr::map_df(sample(1:1000, nperm), run_ranger) %>%
    dplyr::group_by(taxa) %>%
    dplyr::summarise(., average = mean(importance))
  
  model_importance_list <- model_importance %>% dplyr::filter(., average > mean(average)) %>% dplyr::filter(., average > 0) %>% dplyr::pull(., taxa)
  output_sf <- df %>% dplyr::select(., subject_id, feature_of_interest, all_of(model_importance_list))

  return(output_sf)
  
}
