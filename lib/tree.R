## HFE FUNCTIONS

library(janitor, quietly = T, verbose = F, warn.conflicts = F)
library(tidyr, quietly = T, verbose = F, warn.conflicts = F)
library(tibble, quietly = T, verbose = F, warn.conflicts = F)
library(progress, quietly = T, verbose = F, warn.conflicts = F)
library(data.tree, quietly = T, verbose = F, warn.conflicts = F)
library(dplyr, quietly = T, verbose = F, warn.conflicts = F)
options(dplyr.summarise.inform = FALSE)
library(corrr, quietly = T, verbose = F, warn.conflicts = F)
library(tibble, quietly = T, verbose = F, warn.conflicts = F)
library(purrr, quietly = T, verbose = F, warn.conflicts = F)
library(ranger, quietly = T, verbose = F, warn.conflicts = F)
library(vroom, quietly = T, verbose = F, warn.conflicts = F)
library(tidyselect, quietly = T, verbose = F, warn.conflicts = F)
library(recipes, quietly = T, verbose = F, warn.conflicts = F)
library(mikropml, quietly = T, verbose = F, warn.conflicts = F)
suppressPackageStartupMessages(library(ggplot2, quietly = T, verbose = F, warn.conflicts = F))
suppressPackageStartupMessages(library(tidymodels, quietly = T, verbose = F, warn.conflicts = F))
library(fastshap, quietly = T, verbose = F, warn.conflicts = F)
library(shapviz, quietly = T, verbose = F, warn.conflicts = F)
suppressPackageStartupMessages(library(doParallel, quietly = T, verbose = F, warn.conflicts = F))
library(foreach, quietly = T, verbose = F, warn.conflicts = F)

# trim outliers from mean feature abundance calc
# UPDATE: intially we had at 0.02, for an outlier resistant mean
# but if you want a = 0, p = 0 (all features), the trim
# was causing some features to get filtered (they look like
# zero abundance, but really the top 2% of features had non-zero abundances)
trim <- 0.00 

## helper functions ============================================================

## Negate function ("not in"):
`%!in%` <- Negate(`%in%`)

## suppress warnings
options(warn = -1)

## read in metadata  ===========================================================
## rename the subject_identifier to subject_id and
## rename the label to feature_of_interest
## metadata, should be in tab or comma separated format
read_in_metadata <- function(input, subject_identifier, label, feature_type, random_effects, limit_covariates = TRUE, k) {

  cat("\n\n", "Checking for METADATA...", "\n")
  if (file.exists(input) == FALSE) {
    stop("METADATA input not found.")
  }
  cat("\n", paste0("Using ", input, " as METADATA"), "\n")
  
  # read extension to determine file delim
  if (strsplit(basename(input), split = "\\.")[[1]][2] %in% c("tsv","txt")) {
    delim = "\t"
  } else {
    delim = ","
  }
  
  # read in metadata, and select only the subject identifier and
  # feature of interest. Drop NA samples.
  metadata <- suppressMessages(readr::read_delim(file = input, delim = delim)) %>%
    dplyr::rename(., "subject_id" = subject_identifier) %>%
    rename(., "feature_of_interest" = label) %>% 
    janitor::clean_names()
  
  ## make check for NAs and warn user how many rows were dropped
  original_row_count <- nrow(metadata)
  metadata <- metadata %>% tidyr::drop_na()
  new_row_count <- nrow(metadata)
  
  if (original_row_count > new_row_count) {
    warning(paste0((original_row_count - new_row_count), " number of metadata rows were dropped because they contained NAs"), immediate. = TRUE)
    if ((original_row_count - new_row_count) <= 0) {
      stop("All rows were dropped in NA removal.")
    }
  }
  
  ## check and make sure there are not too many metadata columns
  ## we will allow 10 total columns (8 columns of additional covariates)
  if (ncol(metadata) > 10 & limit_covariates) {
    stop("Please only provide subject, label, and up to 8 additional covariates in the metadata file.")
  }
  
  ## notify user what covariates we find, if any
  if (ncol(metadata) > 2) {
    cat("You supplied covariates to consider for the RF competition. They are:\n")
    print(metadata %>% 
            dplyr::select(., -subject_id, -feature_of_interest) %>%
            colnames())
  }
  
  ## this is an effort to clean names for the RF later, which will complain big time
  ## if there are symbols or just numbers in your subject IDs
  metadata$subject_id <- metadata$subject_id %>% janitor::make_clean_names(use_make_names = F)
  
  ## for now, we convert continous response to k levels for a random effects run
  if (random_effects) {
    if (FALSE %in% (c("individual", "time") %in% colnames(metadata))) {
      stop("You specified random effects, you must have metadata columns named individual and time")
    }
    if (feature_type == "numeric") {
      ## this is real ugly code, from SO (https://stackoverflow.com/questions/39906180/consistent-cluster-order-with-kmeans-in-r)
      ## basically it calculates the kmeans and sorts the clusters based on the center means
      ## with those centered means, it applies the cluster index to each sample. Without this,
      ## cluster 2 may repersent the largest values of feature of interest and cluster 3 may repersent the smallest.
      metadata$cluster <- paste0("feature_of_interest_", kmeans(metadata$feature_of_interest, centers = sort(kmeans(metadata$feature_of_interest, centers = k)$centers))$cluster)
      metadata <- metadata %>% dplyr::select(., -feature_of_interest) %>% dplyr::rename(., "feature_of_interest" = "cluster")
    }
  }
  
  return(metadata)
}


## read in hierarchical data =====================================================
## read in data, should be in tab or comma separated format
read_in_hierarchical_data <- function(input, metadata, cores) {
  
  cat("\n", "Checking for DATA...", "\n")
  if (file.exists(input) == FALSE) {
    stop("DATA input not found.")
  }
  cat("\n", paste0("Using ", input, " as DATA"), "\n") 
  
  ## read extension to determine file delim
  if (strsplit(basename(input), split = "\\.")[[1]][2] %in% c("tsv","txt")) {
    delim = "\t"
  } else {
    delim = ","
  }
  
  ## read in hierarchical data, using fast read-in package Vroom.
  ## Useful for large datasets.
  hData <- suppressMessages(vroom::vroom(file = input, delim = delim, skip = 0, 
                                         .name_repair = "minimal", 
                                         num_threads = cores) %>% 
                              dplyr::select(., -any_of(c("NCBI_tax_id", 
                                                         "clade_taxid"))) %>%
                              # clean names so they match metadata and remove
                              # symbols.
                              janitor::clean_names(use_make_names = F))
  
  ## only select columns that are in metadata file, reduce computation if
  ## you have a lot more data than metadata
  
  ## To:do - is this necessary? When we merge to do RF, it will only use the info
  ## that we have metadata for...but abundance and correlation will use all features.
  ## Thought - throwing away samples is ~probably~ not causing a drop in features.
  hData <- hData %>% 
    dplyr::select(., dplyr::any_of(c("clade_name", metadata$subject_id)))
  
  ## try and remove weird symbols in feature names
  hData$clade_name <- stringr::str_replace_all(hData$clade_name, "[^_|[:alnum:]]", "")
  
  ## check and make sure clade_name is the first column in hData
  if ("clade_name" %!in% colnames(hData)) {
    stop("column clade_name not found in input data")
  }
  if (colnames(hData)[1] != "clade_name") {
    hData <- hData %>% dplyr::relocate(., "clade_name")
  }
  
  ## check and make sure there are no NAs in hData
  if (anyNA(hData)) {
    stop("Please remove NAs from your hierarchical data input.")
  }
  
  ## write input to file
  return(as.data.frame(hData))
}

## write separate files to test summary levels =================================
## takes in input from Flatten tree to df function output
# write summarized abundance files for each level except taxa_tree
write_summary_files <- function(input, metadata, output) {
  
  ## write files for all the individual levels
  max_levels <- max(input[["depth"]])
  ## start at 2 to ignore taxa_tree depth (meaningless node)
  levels <- c(1:max_levels)
  
  ## split raw data by pipe symbol into number of expected parts
  count <- 1
  for (i in seq(levels)) {
    if (i == 1) {
      next
    }
    
    ## select different levels 1>i and write them to file
    ## only select features that passed prevalence and abundance thresholds
    file_summary <- input %>%
      dplyr::filter(., depth == i & passed_prevelance == TRUE & passed_abundance == TRUE) %>%
      dplyr::select(., name, 11:dplyr::last_col())
    
    ## there are fringe cases where some levels are named the same but
    ## are part of different clades ie
    ## k_bacteria|p_firmicutes|c_CFG_10299
    ## k_bacteria|p_actinobacteria|c_CFG_10299
    file_summary$name <- file_summary$name %>% janitor::make_clean_names()
    file_summary <- file_summary %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames(., var = "name") %>%
      t() %>%
      as.data.frame() %>%
      tibble::rownames_to_column(., var = "subject_id") %>%
      janitor::clean_names()
    
    ## merge with metadata
    file_merge <- merge(metadata, file_summary, by = "subject_id")
    
    filename <- paste0("_level_", count, ".csv")
    readr::write_delim(x = file_merge, file = paste0(tools::file_path_sans_ext(output), filename), delim = ",")
    count <- count + 1
  }
}

write_taxa_hfe_v1_input_file <- function(input, output) {
  ## write file for TaxaHFE version 1
  version1 <- input %>%
    dplyr::filter(., name != "taxaTree") %>%
    dplyr::select(., pathString, 11:dplyr::last_col()) %>%
    dplyr::rename(., "clade_name" = "pathString") %>%
    tibble::remove_rownames()
  version1$clade_name <- gsub(pattern = "taxaTree\\/", replacement = "", x = version1$clade_name)
  version1$clade_name <- gsub(pattern = "\\/", replacement = "\\|", x = version1$clade_name)
  readr::write_delim(x = version1, file = paste0(tools::file_path_sans_ext(output), "v1_input.csv"), delim = ",")
}

## write files for old_HFE =====================================================
# write old files for the Oudah program
write_oudah_input <- function(input, output) {
  
  ## select "base" (no covariates) metadata file
  metadata <- metadata %>% dplyr::select(., subject_id, feature_of_interest)
  
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
  
  ## get the raw abundance data for the otu.tab input file
  abundance <- input %>%
    dplyr::filter(., depth == max_levels & passed_prevelance == TRUE & passed_abundance == TRUE) %>%
    dplyr::select(., name, 11:dplyr::last_col()) %>%
    tibble::remove_rownames()
  
  ## start writing the taxonomy.tab input file for Oudah HFE  
  taxonomy <- taxonomy[taxonomy[,ncol(taxonomy)] %in% abundance$name, ]
  input_taxa_merge <- merge(taxonomy, abundance, by.x = paste0("L", max_levels), by.y = "name")
  
  input_taxa_merge$index <- (1001:(NROW(input_taxa_merge) + 1000))
  input_taxa_merge$L2 <- "k__Bacteria"
  
  input_taxa_merge <- input_taxa_merge %>%
    dplyr::relocate(., index, dplyr::any_of(c(unlist(paste0("L", c(1:max_levels))))))
  
  ## write OTU to file
  readr::write_delim(x = input_taxa_merge[1:(max_levels)], file = paste0(tools::file_path_sans_ext(output), "_old_hfe_taxa.txt"), col_names = FALSE, delim = "\t")
  readr::write_delim(x = input_taxa_merge %>% dplyr::select(., 1,(max_levels+1):dplyr::last_col()), file = paste0(tools::file_path_sans_ext(output), "_old_hfe_otu.txt"), col_names = FALSE, delim = "\t")
  
  ## write the metadata input for oudah HFE (labels.tab)
  metadata_order <- colnames(input_taxa_merge[,9:NCOL(input_taxa_merge)])
  metadata_list <- metadata %>% dplyr::arrange(match(subject_id, metadata_order)) %>%
    pull(., feature_of_interest)
  ## the first item in that tab seperated list has to be "label
  metadata_list <- as.data.frame(c("label", metadata_list))
  readr::write_delim(x = as.data.frame(t(metadata_list)), file = paste0(tools::file_path_sans_ext(output), "_old_hfe_label.txt"), col_names = FALSE, delim = "\t")
  
}

## get descendant winners =======================================================
# get_descendant_winners goes through the descendants of node, returning a list of all found winners
# if the node level is equal to max_level it can't have any descendant winners by definition, so this function returns an empty list
get_descendant_winners <- function(node, max_level) {
  winners <- list()

  # if the node level is equal to max_level it can't have any descendant winners by definition
  if (node$level == max_level) {
    return(winners)
  }

  for (child in node$children) {
    # if the child is a winner, add the child to list and move to the next child
    if (child$winner) {
      winners <- append(winners, child)
      next
    }

    # otherwise, check the child's children for winners
    winners <- append(winners, get_descendant_winners(child, max_level))
  }

  return(winners)
}

## set the initial values for a leaf ===========================================
# modifies the passed-in node
# node: current leaf node
# row_num: unique row number from original df
# row_vecotr: vector of the numeric data from a single row
# filter_prevalence: passed in filter cutoff for prevalence percentage
# filter_mean_abundance: passed in filter cutoff for mean abundance
initial_leaf_values <- function(node, row_num, row_vector, filter_prevalence, filter_mean_abundance) {
  # a single row dataframe of the abundance data for this row of the df
  node$abundance <- row_vector
  # indicates if the prevalence filter was passed
  node$passed_prevalence_filter <-
    length(node$abundance[node$abundance != 0]) > (length(node$abundance) * filter_prevalence)
  # indicates if the mean abundance filter was passed
  node$passed_mean_abundance_filter <-
    mean(node$abundance, trim = trim) > filter_mean_abundance
  # defaults to be modified later
  node$sf_winner <- FALSE
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
# row_len: length of an abundance vector in the data
# next_row_id: the next id after the last row in the original dataframe
fix_unpopulated_node <- function(node, row_len, next_row_id, filter_prevalence, filter_mean_abundance) {
  # Ignore nodes with abundance or no children
  if (!is.null(node$abundance)) {
    return()
  }

  # sum non-null child abundance vectors
  # if no child abundances are found, the result will be a zero abundance vector
  abundance <- numeric(row_len)
  for (child in node$children) {
    if (is.null(child$abundance)) next
    abundance <- abundance + child$abundance
  }

  # Populate values now that the summed abundance exists
  initial_leaf_values(node, next_row_id, abundance, filter_prevalence, filter_mean_abundance)
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
  pb <- progress::progress_bar$new(format = " Adding nodes to tree [:bar] :percent in :elapsed", total = nrow(df), clear = FALSE, width = 60)

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
    initial_leaf_values(node, row, as.numeric(df[row, 2:ncol(df)]), filter_prevalence, filter_mean_abundance)
  }

  # now that the tree is built, handle unpopulated leaves with the fix_unpopulated_node
  # start the unique id counter at 1 greater than the original df size
  next_row_id <- nrow(df) + 1

  pb2 <- progress::progress_bar$new(format = " Fixing unpopulated nodes [:bar] :percent in :elapsed", total = taxa_tree$totalCount, clear = FALSE, width = 60)

  # traverse the tree and fix the unpopulated nodes
  taxa_tree$Do(function(node) {
    pb2$tick()
    fix_unpopulated_node(node, ncol(df) - 1, next_row_id, filter_prevalence, filter_mean_abundance)
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
compete_node <- function(node, col_names, lowest_level, max_level, corr_threshold, metadata, ncores, feature_type, nperm, random_effects) {
  # skip anything lower than the lowest level (exclusive)
  if (node$level < lowest_level) {
    return()
  }
  # skip anything higher than max_level (exclusive)
  if (node$level > max_level) {
    return()
  }

  ## do not consider children that do not pass abundance and prevalence filters
  if (!node$passed_prevalence_filter || !node$passed_mean_abundance_filter) {
    node$outcomes <- append(node$outcomes, "loss: did not pass filters")
    return()
  }

  # handle no children, this node is the winner
  # also considers a node a winner if it is at the max_level
  if (length(node$children) == 0) {
    node$outcomes <- append(node$outcomes, "win: no children")
    node$winner <- TRUE
    return()
  } else if (node$level == max_level) {
    node$outcomes <- append(node$outcomes, "win: max_level reached")
    node$winner <- TRUE
    return()
  }

  # build dataframe of parent and descendant winners
  # parent is always row 1
  df <- rbind(data.frame(), node$abundance)
  row_names <- c(node$id)

  descendant_winners <- get_descendant_winners(node, max_level)
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
    row_names <- append(row_names, descendant$id)
  }

  rownames(df) <- row_names
  colnames(df) <- col_names

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
  rf_winners <- rf_competition(
    transposed,
    metadata,
    parent_descendent_competition = TRUE,
    "feature_of_interest",
    "subject_id",
    feature_type,
    ncores,
    nperm,
    random_effects
  )

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
compete_all_winners <- function(tree, metadata, col_names, feature_type, nperm, ncores, random_effects) {
  # all vs all competition with winners
  # skipped rows have winner = FALSE so won't appear in this list
  competitors <- get_descendant_winners(tree, tree$height)
  if (length(competitors) == 0) {
    return()
  }

  # all winners into a transposed data frame
  df <- data.frame()
  row_names <- c()
  for (winner in competitors) {
    df <- rbind(df, winner$abundance)
    row_names <- append(row_names, winner$id)
  }

  rownames(df) <- row_names
  colnames(df) <- col_names
  transposed <- as.data.frame(t(df))

  ## return list of winner ids
  rf_winners <- rf_competition(
    transposed,
    metadata,
    parent_descendent_competition = FALSE,
    feature_of_interest = "feature_of_interest",
    subject_identifier = "subject_id",
    feature_type = feature_type,
    ncores = ncores, 
    nperm = nperm,
    random_effects = random_effects
  )


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
      competitor$outcomes <- append(competitor$outcomes, sprintf("win: final rf winner, %s", outcome_str))
      competitor$sf_winner <- TRUE
    } else {
      competitor$outcomes <- append(competitor$outcomes, sprintf("loss: final rf loser, %s", outcome_str))
      competitor$sf_winner <- FALSE
    }
  }
}

## wrapper function to loop through tree nodes =================================
# competes an entries tree
# takes a data.tree root node as input
# modify_tree: determines if the input tree will be modified in place
# col_names: vector of column names from the input df
# lowest_level: lowest level of the tree to consider during operations
#   this level will be compared in a final all-vs-all random forest after the main tree competition
#   defaults to skipping the lowest level (all abundances)
# max_level: determines how deep the descendant competitions will be held
#   defaults to a massive number to allow every descendant
# corr_threshold: the threshold to mark a descendant as highly correlated
# metadata: the metadata associated with the input df that generated the tree
# ncores: the number of cores to use when running the random forest
# disable_super_filter: disables running the final competition
compete_tree <- function(tree, modify_tree = TRUE, col_names, lowest_level = 2, max_level = 1000, corr_threshold, metadata, ncores, feature_type, nperm, disable_super_filter, random_effects) {
  # if not modifying the input tree, create a copy of the tree to perform the competition
  if (!modify_tree) tree <- data.tree::Clone(tree)

  pb <- progress::progress_bar$new(format = " Competing tree [:bar] :percent in :elapsed", total = tree$totalCount, clear = FALSE, width = 60)

  # perform the competition, modifying the tree (which may or may not be a clone of the input)
  tree$Do(
    function(node, col_names, lowest_level, max_level, corr_threshold, metadata, ncores, feature_type, nperm, random_effects) {
      pb$tick()
      compete_node(node, col_names, lowest_level, max_level, corr_threshold, metadata, ncores, feature_type, nperm, random_effects)
    },
    col_names = col_names,
    lowest_level = lowest_level,
    max_level = max_level,
    corr_threshold = corr_threshold,
    metadata = metadata,
    ncores = ncores,
    feature_type = feature_type,
    nperm = nperm,
    traversal = "post-order",
    random_effects = random_effects
  )

  # compete all winners
  # increasing nperm by a factor of 10 to further reduce the variability in the final rf importance scores
  if (disable_super_filter == FALSE) {
    compete_all_winners(
      tree,
      metadata,
      col_names = col_names,
      feature_type = feature_type,
      nperm = nperm * 10,
      ncores = ncores,
      random_effects
    )
  } else {
    cat(" Skipping super filter\n")
  }

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
      dplyr::filter(., .[[2]] >= corrThreshold) %>%
      dplyr::pull(., term)
  )
}

## rf competition function =====================================================
## This is the meat of the ML part of taxaHFE. Input is a dataframe from 
## compete_node() called transposed. This is output from the correlation battle
## of the parent and children who are less correlated than cor_level.
rf_competition <- function(df, metadata, parent_descendent_competition, feature_of_interest = "feature_of_interest", subject_identifier = "subject_id", feature_type, ncores, nperm, random_effects) {
  ## get a list of the covariates in order to remove them from the RF winners
  ## later, so the only RF winners are taxa
  covariates <- metadata %>%
    dplyr::select(., -subject_id, -feature_of_interest) %>%
    colnames()
  # merge node abundance + children abundance with metadata
  merged_data <- merge(df, metadata, by.x = "row.names", by.y = "subject_id")
  merged_data <- merged_data %>% tidyr::drop_na()
  # clean node names so ranger doesnt throw an error
  merged_data <- tibble::column_to_rownames(merged_data, var = "Row.names")
  merged_data <- merged_data %>% janitor::clean_names()
  
  # determine if rf regression or classification should be run
  if (feature_type == "factor") {
    response_formula <- as.formula(paste("as.factor(", feature_of_interest, ") ~ .", sep = ""))
    if (random_effects) {
      merged_data_avg_abund_rf <- prep_re_data(input = merged_data, feature_type = "factor", abund = TRUE)
      merged_data_slope_rf <- prep_re_data(input = merged_data, feature_type = "factor", abund = FALSE)
    }
  } else {
    response_formula <- as.formula(paste("as.numeric(", feature_of_interest, ") ~ .", sep = ""))
    if (random_effects) {
      merged_data_avg_abund_rf <- prep_re_data(input = merged_data, feature_type = "numeric", abund = TRUE)
      merged_data_slope_rf <- prep_re_data(input = merged_data, feature_type = "numeric", abund = FALSE)
    }
  }
  
  # progress bar for the final rf competition
  # will only be incremented/shown if parent_descendent_competition == FALSE
  pb <- progress::progress_bar$new(format = " Competing final winners [:bar] :percent in :elapsed", total = nperm, clear = FALSE, width = 60)
  
  # run ranger, setting parameters such as
  # random seed
  # num.threads number of threads to five ranger
  run_rf <- function(seed, random_effects) {
    if (!parent_descendent_competition) pb$tick()
    
    if (random_effects) {
      ranger_avg_abundance <- ranger::ranger(response_formula, data = merged_data_avg_abund_rf, importance = "impurity_corrected", seed = seed, sample.fraction = 1, replace = TRUE, num.threads = ncores)$variable.importance %>%
        as.data.frame() %>%
        dplyr::rename(., "importance" = ".") %>%
        dplyr::mutate(., importance_rank = rank(importance)) %>%
        dplyr::select(., -importance) %>%
        dplyr::rename(., "importance" = "importance_rank") %>%
        tibble::rownames_to_column(var = "taxa") 
      ranger_slope <- ranger::ranger(response_formula, data = merged_data_slope_rf, importance = "impurity_corrected", seed = seed, sample.fraction = 1, replace = TRUE, num.threads = ncores)$variable.importance %>%
        as.data.frame() %>%
        dplyr::rename(., "importance" = ".") %>%
        dplyr::mutate(., importance_rank = rank(importance)) %>%
        dplyr::select(., -importance) %>%
        dplyr::rename(., "importance" = "importance_rank") %>%
        tibble::rownames_to_column(var = "taxa") 
      ranger_result <- merge(ranger_avg_abundance, ranger_slope, by = "taxa", all = T)
      ranger_result %>% dplyr::mutate(., importance = (importance.x + importance.y)/2) %>% dplyr::select(., -importance.x, -importance.y)
      
    } else {
      ranger::ranger(response_formula, data = merged_data, importance = "impurity_corrected", seed = seed, sample.fraction = 1, replace = TRUE, num.threads = ncores)$variable.importance %>%
        as.data.frame() %>%
        dplyr::rename(., "importance" = ".") %>%
        tibble::rownames_to_column(var = "taxa")
    }
  }
  
  # run the above function across nperm random seeds and average the vip scores
  model_importance <- purrr::map_df(sample(1:1000000, nperm), 
                                    function(seed) run_rf(seed, random_effects = random_effects)) %>%
    dplyr::group_by(taxa) %>%
    dplyr::summarise(., average = mean(importance)) %>%
    dplyr::filter(., taxa %!in% covariates)
  
  # if this is not a parent vs descendent competition
  # return the ids of competitors whose scores meet the following thresholds:
  #   - greater than the average score
  #   - greater than zero
  if (!parent_descendent_competition) {
    return(
      gsub(pattern = "x", replacement = "", x = model_importance %>%
             dplyr::filter(., average > mean(average)) %>%
             dplyr::filter(., average > 0) %>%
             dplyr::pull(., taxa))
    )
  }
  
  # otherwise
  # if top score is the parent, parent wins, else grab the children who
  # beat the parent's score
  # specify the parent column, which is the score to beat
  parentColumn <- janitor::make_clean_names(colnames(df)[1])
  
  ## create a check to make sure the top 2 features are not tied
  tie_check <- model_importance %>% dplyr::arrange(desc(average)) %>% dplyr::slice_head(n = 2)
  if (tie_check$average[1] == tie_check$average[2]) {
    if (parentColumn %in% tie_check$taxa) {
      model_importance$average[model_importance$taxa == parentColumn] <- model_importance$average[model_importance$taxa == parentColumn] + 0.000001
    }
  }
  
  if ((model_importance %>% dplyr::arrange(desc(average)) %>% dplyr::pull(taxa))[1] == parentColumn) {
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
    sf_winner = node$sf_winner,
    rf_loss = node$lost_rf,
    highly_cor = node$highly_correlated,
    passed_prevelance = node$passed_prevalence_filter,
    passed_abundance = node$passed_mean_abundance_filter,
    abundance = data.frame(t(sapply(node$abundance, c))),
    stringsAsFactors = FALSE
  )

  if (length(node$children) > 0) {
    children_df <- do.call(rbind, lapply(node$children, flatten_tree_with_metadata))
    df <- rbind(df, children_df)
  }

  return(df)
}

## massage the output of flattened tree to be less information
prepare_flattened_df <- function(node, metadata, disable_super_filter, col_names) {
  
  ## flatten tree with the supplied metadata
  flattened_df <- flatten_tree_with_metadata(node)
  ## add back in the col names
  colnames(flattened_df)[11:NCOL(flattened_df)] <- col_names
  
  ## only take the features that win or win the superfilter
  if (disable_super_filter) {
    flattened_df <- flattened_df %>%
      dplyr::filter(., winner == TRUE)
  } else {
    flattened_df <- flattened_df %>%
      dplyr::filter(., sf_winner == TRUE)
  }
  
  ## clean pathString names and use these, they will always be unique
  ## downside, longer names
  flattened_df$pathString <- flattened_df$pathString %>% janitor::make_clean_names()
  flattened_df$pathString <- gsub(pattern = "taxa_tree_", replacement = "", x = flattened_df$pathString)
  
  flattened_df <- flattened_df %>%
    dplyr::select(., pathString, 11:dplyr::last_col()) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(., var = "pathString") %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "subject_id")
  
  flattened_df <- merge(metadata, flattened_df, by = "subject_id")
  
  return(flattened_df)
  
}

# write an output file containing the HFE results
write_output_file <- function(flattened_df, metadata, output_location, file_suffix) {
  
  ## clean pathString names and use these, they will always be unique
  ## downside, longer names
  flattened_df$pathString <- flattened_df$pathString %>% janitor::make_clean_names()
  flattened_df$pathString <- gsub(pattern = "taxa_tree_", replacement = "", x = flattened_df$pathString)
  
  output <- flattened_df %>%
    dplyr::select(., pathString, 11:dplyr::last_col()) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(., var = "pathString") %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "subject_id")

  output <- merge(metadata, output, by = "subject_id")
  readr::write_delim(file = paste0(tools::file_path_sans_ext(output_location), file_suffix), x = output, delim = ",")
}

# generate the outputs
# if disable_super_filter is TRUE, the super filter competition wasn't run, and only one output will be generated
# if both_outputs is FALSE, one file will be produced with the final level of competition that occurred
generate_outputs <- function(tree, metadata, col_names, output_location, disable_super_filter, write_both_outputs, write_old_files, write_flattened_df_backup, ncores) {
  # flatten tree back into metadata, and assign original column names from hData to the sample columns
  flattened_df <- flatten_tree_with_metadata(tree)
  colnames(flattened_df)[11:NCOL(flattened_df)] <- col_names
  
  ## filter to only winners and clean names in case of duplicate
  ## also further filtering for sf winners
  flattened_winners <- flattened_df %>%
    dplyr::filter(., winner == TRUE)

  flattened_sf_winners <- flattened_winners %>%
    dplyr::filter(., sf_winner == TRUE)

  ## if super filter is disabled, write the flattened_winners as standard output
  ## otherwise write the super filter winners as standard output
  if (disable_super_filter == TRUE) {
    write_output_file(flattened_winners, metadata, output_location, ".csv")
  } else {
    write_output_file(flattened_sf_winners, metadata, output_location, ".csv")
  }

  ## also write the non-sf output if both outputs are requested
  if (write_both_outputs == TRUE && disable_super_filter == FALSE) {
    write_output_file(flattened_winners, metadata, output_location, "_no_sf.csv")
  }

  cat(" Features (no super filter): ", nrow(flattened_winners), "\n")
  if (disable_super_filter != TRUE) {
    cat("\n Features (super filter): ", nrow(flattened_sf_winners), "\n")
  }

  ## write old files
  if (write_old_files == TRUE) {
    cat("\n", "###########################\n", "Writing old files...\n", "###########################\n\n")

    write_summary_files(input = flattened_df, metadata = metadata, output = output_location)
    write_oudah_input(input = flattened_df, output = output_location)
  }

  ## save flattened DF to come back to
  if (write_flattened_df_backup == TRUE) {
    vroom::vroom_write(
      x = flattened_df,
      file = paste0(tools::file_path_sans_ext(output_location), "_raw_data.tsv.gz"),
      num_threads = ncores
    )
  }
}

## simple function to store objects into diet_ml_inputs list, with
## custom attributes that keep track of
## 1. program_method (taxa_hfe_ml, summarized levels, etc.)
## 2. whether superfilter was used
## 3. is it a train or test object
## 4. what summarized level is it
## 5. what random seed was run
store_diet_ml_inputs <- function(target_list, object, super_filter, method, train_test_attr, level_n, seed) {
  
  target_name <- paste0(method, "_", super_filter, "_", train_test_attr, "_", level_n)
  ## add target object to target_list with object_name
  target_list[target_name] <- list(object)
  ## add attribute that tells us what program it came from
  attr(target_list[[target_name]], "program_method") <- method
  ## add attribute that tells us if superfilter was used
  attr(target_list[[target_name]], "superfilter") <- super_filter
  ## add attribute about training or testing
  attr(target_list[[target_name]], "train_test_attr") <- train_test_attr
  ## add attribute about what level the data was summarized
  attr(target_list[[target_name]], "level") <- level_n
  ## add attribute about seed was used
  attr(target_list[[target_name]], "seed") <- seed
  
  return(target_list)
}

## same function as write_summary_files() except it just 
## adds these objects to diet_ml_inputs list.
generate_summary_files <- function(input, metadata, target_list, object, 
                                   disable_super_filter, seed) {
  
  ## write files for all the individual levels
  max_levels <- max(input[["depth"]])
  ## start at 2 to ignore taxa_tree depth (meaningless node)
  levels <- c(1:max_levels)
  
  ## split raw data by pipe symbol into number of expected parts
  count <- 1
  for (i in seq(levels)) {
    if (i == 1) {
      next
    }
    
    ## select different levels 1>i and write them to file
    ## only select features that passed prevalence and abundance thresholds
    flattened_df <- input %>% dplyr::filter(., depth == i & passed_prevelance == TRUE & passed_abundance == TRUE) 
    
    ## clean pathString names and use these, they will always be unique
    ## downside, longer names
    flattened_df$pathString <- flattened_df$pathString %>% janitor::make_clean_names()
    flattened_df$pathString <- gsub(pattern = "taxa_tree_", replacement = "", x = flattened_df$pathString)
    
    output <- flattened_df %>%
      dplyr::select(., pathString, 11:dplyr::last_col()) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames(., var = "pathString") %>%
      t() %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "subject_id")
    
    ## merge with metadata
    level <- merge(metadata, output, by = "subject_id")
    
    diet_ml_inputs <- store_diet_ml_inputs(target_list = diet_ml_inputs,
                                          object = level,
                                          super_filter = NA,
                                          method = "summarized_level",
                                          train_test_attr = NA,
                                          level_n = count,
                                          seed = seed
    )
    
    count <- count + 1
  }
  
  return(diet_ml_inputs)
}

## loops over dietML_inputs list and checks if there is a train, test
## version of the object (df). If there isnt, creates 2 new objects in the
## list that is a test and a train object of that original object. The split
## is informed from the tr_te_split code in the run_file.R
split_train_data <- function(target_list, attribute_name, seed) {
  # Initialize an empty vector to store indices with NA values
  na_indices <- integer(0)
  
  # Loop through the list
  for (i in seq_along(target_list)) {
    # Get the current item
    item <- target_list[[i]]
    
    # Check if the item has the attribute
    if (!is.null(attr(item, attribute_name))) {
      # Get the value of the attribute
      attr_value <- attr(item, attribute_name)
      
      # Check if the attribute value is NA
      if (is.na(attr_value)) {
        # Append the index to the na_indices vector
        na_indices <- c(na_indices, i)
      }
      }
    }
  
  for (missing_train_index in na_indices) {
    temp_train <- diet_ml_inputs[[missing_train_index]] %>% as.data.frame() %>% dplyr::filter(., subject_id %in% train_metadata$subject_id)
    diet_ml_inputs <- store_diet_ml_inputs(target_list = diet_ml_inputs,
                                          object = temp_train,
                                          super_filter = attributes(diet_ml_inputs[[missing_train_index]])$superfilter,
                                          method = attributes(diet_ml_inputs[[missing_train_index]])$program_method,
                                          train_test_attr = "train",
                                          level_n = attributes(diet_ml_inputs[[missing_train_index]])$level,
                                          seed = seed
    )
    
    temp_test <- diet_ml_inputs[[missing_train_index]] %>% as.data.frame() %>% dplyr::filter(., subject_id %in% test_metadata$subject_id)
    diet_ml_inputs <- store_diet_ml_inputs(target_list = diet_ml_inputs,
                                          object = temp_test,
                                          super_filter = attributes(diet_ml_inputs[[missing_train_index]])$superfilter,
                                          method = attributes(diet_ml_inputs[[missing_train_index]])$program_method,
                                          train_test_attr = "test",
                                          level_n = attributes(diet_ml_inputs[[missing_train_index]])$level,
                                          seed = seed
    )
  }
  
  return(diet_ml_inputs)
  
}

## writes every object in dietML_inputs to file
write_list_to_csv <- function(target_list, directory = ".") {
  # Check if the provided directory exists
  if (!dir.exists(directory)) {
    stop("The specified directory does not exist.")
  }
  
  # Get the names of the list items
  names_lst <- names(target_list)
  
  # Check if the list has names
  if (is.null(names_lst)) {
    stop("The list must have names for the objects.")
  }
  
  # Loop through the list and write each object to a CSV file
  for (i in seq_along(target_list)) {
    # Get the current item and its name
    item <- target_list[[i]]
    item_name <- names_lst[i]
    
    # Check for the "train_test_attr" attribute and if it is non-NA and non-null
    train_test_attr <- attr(item, "train_test_attr")
    if (is.null(train_test_attr) || is.na(train_test_attr)) {
      message(paste("Skipping", item_name, ": 'train_test_attr' is NULL or NA."))
      next
    }
    
    # Create a filename for the current item
    filename <- file.path(directory, paste0(item_name, ".csv"))
    
    # Check if the item is a data.frame or matrix
    if (is.data.frame(item) || is.matrix(item)) {
      # Write the item to a CSV file
      readr::write_csv(item, filename)
    } else {
      # Print a warning if the item is not a data.frame or matrix
      warning(paste("Item", item_name, "is not a data.frame or matrix and was not written to a CSV file."))
    }
  }
  
  # Inform the user that the process is complete
  message("Finished writing objects to CSV files.")
}

## create a dataframe of attributes I care about (attr_to_return list) from the
## diet_ml_inputs list, which i can use to pass to dietML
extract_attributes <- function(items_list) {
  # Initialize an empty list to store attributes for each item
  attr_list <- list()
  ## create a list of attributes we want to pull into a dataframe
  attr_to_return <- c("program_method", "superfilter", "train_test_attr", "level", "seed")
  
  # Loop through each item in the list
  for (i in seq_along(items_list)) {
    # Get attributes of the current item
    item_attr <- attributes(items_list[[i]])
    item_attr <- subset(item_attr, names(item_attr) %in% attr_to_return)
    item_attr <- append(item_attr, values = c("name" =  names(diet_ml_inputs[i])))
    
    # Add the item number or name for reference
    if (is.null(item_attr)) {
      item_attr <- list(item_name = paste0("Item_", i))
    } else {
      item_attr$item_name <- paste0("Item_", i)
    }
    
    # Store the attributes in the list
    attr_list[[i]] <- as.data.frame(item_attr, stringsAsFactors = FALSE)
  }
  
  # Combine all attributes into one dataframe
  combined_df <- do.call(rbind, attr_list)
  # remove any objects that do not have a train_test_attr. These are the 
  # objects that existed in the list prior to running the split_train_data() 
  # funtion.
  combined_df <- combined_df %>% dplyr::filter(., train_test_attr != "")
  # create a general name for each method. ie instead of summarized_level_NA_train_1,
  # the general name should be summarized_level_1
  combined_df$general_name <- gsub(pattern = "_train|_test|_NA",replacement = "", x = combined_df$name)
  
  return(combined_df)
}

## run dietML based on diet_ml_input_df
run_diet_ml <- function(input_df, n_repeat, feature_type, seed, train, test, model, program, random_effects, folds, cor_level, ncores, tune_length, tune_stop, tune_time, metric, label, output, shap) {
  
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
                     program = program, 
                     model = model, 
                     seed = seed, 
                     random_effects, 
                     folds = folds, 
                     cor_level = cor_level, 
                     ncores = ncores, 
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

## prep data for random effects random forest
## input needs to be formated as merged data prior to RF
## if abund = TRUE, calc avg abund, else calc avg slope
prep_re_data <- function(input, feature_type, abund) {
  
  ## onehot encode any factors (only factors should be those from covariates)
  merged_data_hotencode <- input %>%
    recipes::recipe(~ .) %>% 
    recipes::step_dummy(recipes::all_nominal_predictors(), -feature_of_interest, -individual, -time) %>% 
    recipes::prep() %>% 
    recipes::bake(input) 
  
  
  ## if response is factor, average abundance across label levels
  ## also average slope of features within individual and label levels
  if (feature_type == "factor") {
    ## abund analysis ##
    
    ## average abundance within level and feature_of_interest
    merged_data_avg_abund <- merged_data_hotencode %>% 
      dplyr::group_by(., individual, feature_of_interest) %>% dplyr::summarise(., across(where(is.numeric), mean)) %>%
      dplyr::ungroup() 
    
    ## remove individual and time columns from dataset for RF
    merged_data_avg_abund_rf <- merged_data_avg_abund %>% dplyr::select(., -individual, -time)
    
    if (abund) {
      return(merged_data_avg_abund_rf)
    } else {
      ## slope analysis ##
      
      ## all numerical data
      numerical_cols <- input %>% select(where(is.numeric)) %>% colnames()
      
      ## get the separate averages for the one hot encoded factors from above
      one_hot_encoded_cov <- merged_data_avg_abund %>%
        dplyr::select(., -dplyr::any_of(numerical_cols), feature_of_interest, individual)
      
      ## get the slopes of the features within feature of interest, across time
      merged_data_slope <- merged_data_hotencode %>% 
        dplyr::select(., dplyr::any_of(numerical_cols), individual, feature_of_interest) %>%
        group_by(individual, feature_of_interest) %>%
        summarise(across(is.numeric,
                         list(slope = ~lm(. ~ as.numeric(time))$coef[2]))) %>%
        dplyr::mutate(across(everything(), ~replace_na(.x, 0))) %>%
        dplyr::select(., -time_slope)
      
      ## remove slope from colname
      colnames(merged_data_slope) <- gsub(pattern = "_slope", replacement = "", x = colnames(merged_data_slope))
      
      ## merge back with factor covariates
      merged_data_slope_rf <- merge(one_hot_encoded_cov, merged_data_slope, by = c("individual", "feature_of_interest"), all = TRUE)
      merged_data_slope_rf <- merged_data_slope_rf %>% dplyr::select(., -dplyr::any_of(c("individual", "time")))
      
      return(merged_data_slope_rf)
    }
  } 
}
  
pass_to_dietML <- function(train, test, model, program, seed, random_effects, folds, cor_level, ncores, tune_length, tune_stop, tune_time, metric, label, output, feature_type, shap) {
  
  ## check and make sure diet_ml_input_df exists
  if (!exists("diet_ml_input_df") | nrow(diet_ml_input_df) < 1) {
    stop(paste0("Nothing passed to dietML."))
  } 
  
  ## check for outdir and make if not there
  if (dir.exists(paste0(dirname(opts$OUTPUT), "/ml_analysis")) != TRUE) {
    dir.create(path = paste0(dirname(opts$OUTPUT), "/ml_analysis"))
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
  results_df <- run_null_model(train = train, test = test, seed = seed, type = type, cor_level = cor_level, random_effects = random_effects, output = output)
  
  ## if specified, run random forest
  if (model == "rf") {
    shap_inputs <- run_dietML_ranger(
      train = train,
      test = test,
      seed = seed,
      random_effects = random_effects,
      folds = folds,
      cor_level = cor_level,
      ncores = ncores,
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
  
  ## if specified, run elastic net
  if (model == "enet") {
    shap_inputs <- run_dietML_enet(
      train = train,
      test = test,
      seed = seed,
      random_effects = random_effects,
      folds = folds,
      cor_level = cor_level,
      ncores = ncores,
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
                  ncores = ncores)
  }
}

run_dietML_ranger <- function(train, test, seed, random_effects, folds, cor_level, ncores, tune_length, tune_stop, tune_time, metric, label, model, program, output, type, null_results) {
  
  ## check and make sure results_df has been created from run_null_model,
  ## needed for this function downstream
  if (!exists("null_results")) {
    stop("Null model was not run. Cannot complete model evaluation.")
  }
  
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
  folds <- rsample::vfold_cv(train, v = as.numeric(folds), strata = feature_of_interest, repeats = 3)
  
  ## recipe
  
  ## specify recipe (this is like the pre-process work)
  diet_ml_recipe <- 
    recipes::recipe(feature_of_interest ~ ., data = train) %>% 
    recipes::update_role(tidyr::any_of("subject_id"), new_role = "ID") %>% 
    recipes::step_dummy(recipes::all_nominal_predictors()) %>%
    recipes::step_corr(recipes::all_numeric_predictors(), threshold = as.numeric(cor_level), use = "everything") %>%
    recipes::step_zv(recipes::all_predictors())
  
  ## ML engine
  
  ## specify ML model and engine 
  initial_mod <- parsnip::rand_forest(mode = type, 
                                      mtry = tune(),
                                      trees = 1500,
                                      min_n = tune()) %>%
    parsnip::set_engine("ranger", 
                        num.threads = 1,
                        importance = "none")
  
  initial_mod %>% parsnip::translate()
  
  ## workflow
  
  ## define workflow
  diet_ml_wflow <- 
    workflows::workflow() %>% 
    workflows::add_model(initial_mod) %>% 
    workflows::add_recipe(diet_ml_recipe)  
  
  ## set up parallel jobs
  ## remove any doParallel job setups that may have
  ## unneccessarily hung around
  unregister_dopar()
  
  ## register parallel cluster
  cl <- parallel::makePSOCKcluster(as.numeric(ncores))
  doParallel::registerDoParallel(cl)
  
  ## hyperparameters
  
  ## define the hyper parameter set
  diet_ml_param_set <- parsnip::extract_parameter_set_dials(diet_ml_wflow)
  
  ## for random forests, set mtry to max features after correlation
  ## co-correlate features at specified threshold (get upper limit of mtry)
  training_cor <- mikropml:::group_correlated_features(train %>% dplyr::select(where(is.numeric)) %>% dplyr::select(., -dplyr::any_of(c("feature_of_interest", "subject_id"))), 
                                                       corr_thresh = as.numeric(cor_level), group_neg_corr = T)
  
  ## make dataframe of what is correlated at specified threshold.
  training_cor <- as.data.frame(training_cor) %>% 
    tidyr::separate(., col = training_cor, into = c("keep", "co_correlated"), sep = "\\|", extra = "merge")
  
  ## set mtry to max features after correlation
  diet_ml_param_set <- 
    diet_ml_param_set %>% 
    # Pick an upper bound for mtry: 
    recipes::update(mtry = mtry(range(c(2, round((NROW(training_cor) * 0.9), digits = 0)))), 
                    min_n = min_n(range(c(2, nrow(test)))))
  
  ## set up hyper parameter search
  if (type == "classification") {
    
    search_res <-
      diet_ml_wflow %>% 
      tune::tune_bayes(
        resamples = folds,
        # To use non-default parameter ranges
        param_info = diet_ml_param_set,
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
                                      seed = as.numeric(seed))
      )
    
  } else if (type == "regression") {
    
    search_res <-
      diet_ml_wflow %>% 
      tune::tune_bayes(
        resamples = folds,
        # To use non-default parameter ranges
        param_info = diet_ml_param_set,
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
                                      seed = as.numeric(seed))
      )
  }
  
  tune::show_best(search_res)
  
  ## stop parallel jobs
  parallel::stopCluster(cl)
  ## remove any doParallel job setups that may have
  ## unneccessarily hung around
  unregister_dopar()
  
  ## fit best model
  
  ## get the best parameters from tuning
  best_mod <- 
    search_res %>% 
    tune::select_best(metric = metric)
  
  ## create the last model based on best parameters
  last_best_mod <- 
    parsnip::rand_forest(mtry = best_mod$mtry, min_n = best_mod$min_n) %>% 
    parsnip::set_engine("ranger", num.threads = as.numeric(ncores), importance = "none") %>% 
    parsnip::set_mode(type)
  
  ## update workflow with best model
  best_tidy_workflow <- 
    diet_ml_wflow %>% 
    workflows::update_model(last_best_mod)
  
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
  readr::write_csv(x = full_results, file = paste0(dirname(output), "/ml_analysis/ml_results.csv"), 
                   append = T, col_names = !file.exists(paste0(dirname(output), "/ml_analysis/ml_results.csv")))
  
  ## load up list for shap analysis
  shap_inputs <- list("split_from_data_frame" = split_from_data_frame, "diet_ml_recipe" = diet_ml_recipe, "best_tidy_workflow" = best_tidy_workflow)
  return(shap_inputs)
  
}

run_null_model <- function(train, test, seed, type, cor_level, random_effects, output) {
  
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
  diet_ml_recipe <- 
    recipes::recipe(feature_of_interest ~ ., data = train) %>% 
    recipes::update_role(tidyr::any_of("subject_id"), new_role = "ID") %>%
    recipes::step_dummy(recipes::all_nominal_predictors()) %>%
    recipes::step_corr(all_numeric_predictors(), threshold = cor_level) %>%
    recipes::step_zv(all_predictors())
  
  
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
    workflows::add_recipe(diet_ml_recipe)  
  
  
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
  readr::write_csv(x = results_df, file =paste0(dirname(output), "/ml_analysis/dummy_model_results.csv"), 
                   append = T, col_names = !file.exists(paste0(dirname(output), "/ml_analysis/dummy_model_results.csv")))
  
  ## return results_df because that is what the other models need (ranger, enet)
  return(results_df)
  
}

shap_analysis <- function(label, output, model, filename, shap_inputs, train, test, type, ncores) {
  
  # --- Load SHAP inputs ---
  split_from_data_frame <- shap_inputs$split_from_data_frame
  best_tidy_workflow <- shap_inputs$best_tidy_workflow
  diet_ml_recipe <- shap_inputs$diet_ml_recipe
  
  # --- Setup ---
  shap_plot_env <- new.env()
  shap.error.occured <- FALSE
  error_message <- NULL
  output_dir <- file.path(dirname(output), "ml_analysis")
  
  # --- Define prediction wrapper (pfun) ---
  pfun <- NULL
  if (length(levels(as.factor(split_from_data_frame$data$feature_of_interest))) == 2) {
    # Binary classification
    pfun <- function(object, newdata) {
      preds <- predict(object, data = newdata)$predictions
      class_level <- levels(as.factor(split_from_data_frame$data$feature_of_interest))[1]
      #print("Binary classification prediction:") ## print commands in pfun() there for debugging
      #print(head(preds))
      #print(paste("Using class level:", class_level))
      return(preds[, class_level])
    }
  } else if (type == "regression" && model == "rf") {
    # Regression
    pfun <- function(object, newdata) {
      preds <- predict(object, data = newdata)$predictions
      #print("Regression prediction:")
      #print(head(preds))
      return(preds)
    }
  }
  
  if (is.null(pfun)) {
    message("Error: Could not define prediction function (pfun). Check model and type inputs.")
    shap.error.occured <- TRUE
  } else {
    # --- SHAP analysis block ---
    result <- tryCatch({
      
      shap_data_subsets <- list(list(split_from_data_frame$data, "full"), list(train, "train"), list(test, "test"))
      
      if (ncores > 1) {
        parallel_shap <- TRUE
        doParallel::registerDoParallel(cores = ncores)
      } else {
        parallel_shap <- FALSE
      }
      
      for (i in seq_along(shap_data_subsets)) {
        # Fit the model
        best_workflow <- parsnip::fit(best_tidy_workflow, shap_data_subsets[[i]][[1]])
        best_workflow_mod <- workflows::extract_fit_parsnip(best_workflow)
        
        # Prepare data
        shap_data <- recipes::prep(diet_ml_recipe, shap_data_subsets[[i]][[1]]) %>%
          recipes::juice() %>%
          dplyr::select(-feature_of_interest, -subject_id)
        
        # Compute SHAP values
        shap_explanations <- fastshap::explain(
          object = best_workflow_mod$fit,
          X = shap_data,
          pred_wrapper = pfun,
          nsim = 500,
          adjust = TRUE,
          parallel = parallel_shap
        )
        
        # SHAP object for plotting
        sv <- shapviz::shapviz(shap_explanations, X = shap_data)
        
        # Generate and save plot
        plot <- shap_plot(
          sv = sv,
          label = label,
          data_subset_label = shap_data_subsets[[i]][[2]],
          split_from_data_frame = split_from_data_frame,
          filename = filename,
          output_dir = output_dir,
          data_subset_index = i
        )
        
        ## these next lines save objects so that shapviz and be re-plotted
        assign(paste0("shap_data_", shap_data_subsets[[i]][[2]]), shap_data, envir = shap_plot_env)
        assign(paste0("shap_explanations_", shap_data_subsets[[i]][[2]]), shap_explanations, envir = shap_plot_env)
        assign(paste0("sv_", shap_data_subsets[[i]][[2]]), sv, envir = shap_plot_env)
        assign(paste0("plot_", shap_data_subsets[[i]][[2]]), plot, envir = shap_plot_env)
        
      }
      
      ## more things to save, outside the for-loop
      assign("split_from_data_frame", split_from_data_frame, envir = shap_plot_env)
      assign("label", label, envir = shap_plot_env)
      
    }, error = function(e) {
      shap.error.occured <<- TRUE
      error_message <<- e$message
      NULL
    })
  }
  
  # --- Save and return results ---
  if (shap.error.occured) {
    message("❌ SHAP analysis could not be completed.")
    if (!is.null(error_message)) {
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
    data_subset_index
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
    labs(x = paste0(
      "predictive of ",
      class_levels[2],
      " < SHAP > predictive of ",
      class_levels[1]
    )) +
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
## get used in this codebase, but its too got not to exist somewhere. 
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

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
