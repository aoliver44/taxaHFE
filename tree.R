## HFE FUNCTIONS
## v2.11

library(janitor, quietly = T, verbose = F, warn.conflicts = F)
library(tidyr, quietly = T, verbose = F, warn.conflicts = F)
library(tibble, quietly = T, verbose = F, warn.conflicts = F)
library(progress, quietly = T, verbose = F, warn.conflicts = F)
library(data.tree, quietly = T, verbose = F, warn.conflicts = F)
library(dplyr, quietly = T, verbose = F, warn.conflicts = F)
library(corrr, quietly = T, verbose = F, warn.conflicts = F)
library(tibble, quietly = T, verbose = F, warn.conflicts = F)
library(purrr, quietly = T, verbose = F, warn.conflicts = F)
library(ranger, quietly = T, verbose = F, warn.conflicts = F)
library(vroom, quietly = T, verbose = F, warn.conflicts = F)
library(tidyselect, quietly = T, verbose = F, warn.conflicts = F)


## set random seed, defaults to system time
set_seed_func <- function(seed) {
  if (!is.null(seed)) {
    opt$seed <<- as.numeric(opt$seed)
    set.seed(seed)
  } else {
    message('No random seed set. Using system time!')
    opt$seed <<- as.numeric(Sys.time())
    set.seed(opt$seed)
  }
}

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
read_in_metadata <- function(input, subject_identifier, label) {
  
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
    dplyr::filter(., !is.na(subject_id)) %>%
    janitor::clean_names()
  
  ## check and make sure there are not too many metadata columns
  ## we will allow 10 total columns (8 columns of additional covariates)
  if (ncol(metadata) > 10) {
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

write_taxaHFEv1_input_file <- function(input, output) {
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

## get descenant winners =======================================================
# get_descendant_winners goes through the descendants of node, returning a list of all found winners
# maxDepth defines how deep the winner function will go to find a winner
get_descendant_winners <- function(node, max_level) {
  winners <- list()

  # if maxDepth is zero, this is the bottom
  # return an empty list as no further generations will be considered
  if (max_level == 0) {
    return(winners)
  }

  for (child in node$children) {
    # if the child is a winner, add the child to list and move to the next child
    if (child$winner) {
      winners <- append(winners, child)
      next
    }

    # otherwise, check the child's children for winners
    winners <- append(winners, get_descendant_winners(child, max_level - 1))
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
compete_node <- function(node, col_names, lowest_level, max_level, corr_threshold, metadata, ncores, feature_type, nperm) {
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
    nperm
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
compete_all_winners <- function(tree, metadata, col_names, feature_type, nperm, ncores) {
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
    ncores = ncores, nperm = nperm
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
compete_tree <- function(tree, modify_tree = TRUE, col_names, lowest_level = 2, max_level = 1000, corr_threshold, metadata, ncores, feature_type, nperm, disable_super_filter) {
  # if not modifying the input tree, create a copy of the tree to perform the competition
  if (!modify_tree) tree <- data.tree::Clone(tree)

  pb <- progress::progress_bar$new(format = " Competing tree [:bar] :percent in :elapsed", total = tree$totalCount, clear = FALSE, width = 60)

  # perform the competition, modifying the tree (which may or may not be a clone of the input)
  tree$Do(
    function(node, col_names, lowest_level, max_level, corr_threshold, metadata, ncores, feature_type, nperm) {
      pb$tick()
      compete_node(node, col_names, lowest_level, max_level, corr_threshold, metadata, ncores, feature_type, nperm)
    },
    col_names = col_names,
    lowest_level = lowest_level,
    max_level = max_level,
    corr_threshold = corr_threshold,
    metadata = metadata,
    ncores = ncores,
    feature_type = feature_type,
    nperm = nperm,
    traversal = "post-order"
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
      ncores = ncores
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
# TODO: document these inputs
rf_competition <- function(df, metadata, parent_descendent_competition, feature_of_interest = "feature_of_interest", subject_identifier = "subject_id", feature_type, ncores, nperm) {
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
  data_colnames <- colnames(merged_data)
  merged_data <- merged_data %>% janitor::clean_names()

  # determine if rf regression or classification should be run
  if (feature_type == "factor") {
    response_formula <- as.formula(paste("as.factor(", feature_of_interest, ") ~ .", sep = ""))
  } else {
    response_formula <- as.formula(paste("as.numeric(", feature_of_interest, ") ~ .", sep = ""))
  }

  # progress bar for the final rf competition
  # will only be incremented/shown if parent_descendent_competition == FALSE
  pb <- progress::progress_bar$new(format = " Competing final winners [:bar] :percent in :elapsed", total = nperm, clear = FALSE, width = 60)

  # run ranger, setting parameters such as
  # random seed
  # num.threads number of threads to five ranger
  run_ranger <- function(seed) {
    if (!parent_descendent_competition) pb$tick()

    ranger::ranger(response_formula, data = merged_data, importance = "impurity_corrected", seed = seed, sample.fraction = 1, replace = TRUE, num.threads = ncores)$variable.importance %>%
      as.data.frame() %>%
      dplyr::rename(., "importance" = ".") %>%
      tibble::rownames_to_column(var = "taxa")
  }

  # run the above function across nperm random seeds and average the vip scores
  model_importance <- purrr::map_df(sample(1:1000000, nperm), run_ranger) %>%
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

  ## clean up the names in case there are weird issues of duplicates
  flattened_df$name <- janitor::make_clean_names(flattened_df$name)
  
  ## get rid of all the run information
  flattened_df <- flattened_df %>%
    dplyr::select(., name, 11:dplyr::last_col()) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(., var = "name") %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "subject_id")
  
  flattened_df <- merge(metadata, flattened_df, by = "subject_id")
  
  return(flattened_df)
  
}

# write an output file containing the HFE results
write_output_file <- function(flattened_df, metadata, output_location, file_suffix) {
  output <- flattened_df %>%
    dplyr::select(., name, 11:dplyr::last_col()) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(., var = "name") %>%
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

  ## write old files  ============================================================
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

## simple function to store objects into dietML_inputs list, with
## custom attributes that keep track of
## 1. program_method (taxaHFE_ML, summarized levels, etc.)
## 2. whether superfilter was used
## 3. is it a train or test object
## 4. what summarized level is it
## 5. what random seed was run
store_dietML_inputs <- function(target_list, object, super_filter, method, train_test_attr, level_n, seed) {
  
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
## adds these objects to dietML_inputs list.
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
    level <- merge(metadata, file_summary, by = "subject_id")
    
    dietML_inputs <<- store_dietML_inputs(target_list = dietML_inputs,
                                          object = level,
                                          super_filter = NA,
                                          method = "summarized_level",
                                          train_test_attr = NA,
                                          level_n = count,
                                          seed = seed
    )
    
    count <- count + 1
  }
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
    temp_train <- dietML_inputs[[missing_train_index]] %>% as.data.frame() %>% dplyr::filter(., subject_id %in% train_metadata$subject_id)
    dietML_inputs <<- store_dietML_inputs(target_list = dietML_inputs,
                                          object = temp_train,
                                          super_filter = attributes(dietML_inputs[[missing_train_index]])$superfilter,
                                          method = attributes(dietML_inputs[[missing_train_index]])$program_method,
                                          train_test_attr = "train",
                                          level_n = attributes(dietML_inputs[[missing_train_index]])$level,
                                          seed = seed
    )
    
    temp_test <- dietML_inputs[[missing_train_index]] %>% as.data.frame() %>% dplyr::filter(., subject_id %in% test_metadata$subject_id)
    dietML_inputs <<- store_dietML_inputs(target_list = dietML_inputs,
                                          object = temp_test,
                                          super_filter = attributes(dietML_inputs[[missing_train_index]])$superfilter,
                                          method = attributes(dietML_inputs[[missing_train_index]])$program_method,
                                          train_test_attr = "test",
                                          level_n = attributes(dietML_inputs[[missing_train_index]])$level,
                                          seed = seed
    )
  }
  
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
## dietML_inputs list, which i can use to pass to dietML
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
    item_attr <- append(item_attr, values = c("name" =  names(dietML_inputs[i])))
    
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

## run dietML based on dietML_input_df
run_dietML <- function(input_df, n_repeat) {
  
  for (seed in sample(1:100000000, replace = F, size = opt$permute)) {
    if (opt$permute > 1) {
      ## same idea as below, need to change the random seed and it needs
      ## to be a persistant change through sourcing the dietML script.
      ## again, we can create a list of seeds that can be passed through
      opt$seed <<- seed
    }
    for (dML_input in unique(input_df[["general_name"]])) {
      
      ## I dont know how else to create these objects and source script within
      ## this function, without assigning to global env
      ## DietML.R is looking for train_data and test_data in global env. I suppose
      ## we could create a list of these objects and tweak dietML to take in
      ## the list.
      
      train_data <<- dietML_inputs[[input_df %>% 
                                     dplyr::filter(., general_name == dML_input & train_test_attr == "train") %>% 
                                     dplyr::pull(name)]]
      
      test_data <<- dietML_inputs[[input_df %>% 
                                    dplyr::filter(., general_name == dML_input & train_test_attr == "test") %>% 
                                    dplyr::pull(name)]]
      
      ## keep track of what method is being passed to dietML
      ## this gets printed in the results file
      opt$program <<- dML_input
      #TODO: make this a function that takes instead of a source. Make "program" and 
      #"seed" and whatnot into function arguments and not something that gets written
      #to opt!!! BANISH UNNECESSARY GLOBAL VARS!!!
      source("/home/docker/taxaHFE-ML/dietML.R")
    }
  }
}
