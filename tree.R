library(data.tree)

hData <- read.delim("/home/docker/example_inputs/microbiome_data.txt")

hTree <- Node$new("taxaTree", id=0)

filterPrevalence <- 0.01
filterMeanAbundance <- 0.0001
corrThreshold <- 0.95
trim <- 0.02
options(width=150)

# go through the descendants of node, returning a list of all found winners
# maxDepth defines how deep the winner function will go to find a winner
getDescendantWinners <- function(node, maxDepth) {
    # if no children or maxDepth is zero, return empty list
    if (length(node$children) == 0 | maxDepth == 0) {
        return(list())
    }

    winnerList <- list()
    for (i in 1:length(node$children)) {
        if (node$children[[i]]$winner) {
            winnerList <- append(winnerList, node$children[[i]])
            next
        }
        winnerList <- append(winnerList, getDescendantWinners(node$children[[i]], maxDepth-1))
    }

    return(winnerList)
}

calculateCorrelation <- function(df, corrThreshold) correlated <- {
    parentColumn <- colnames(df)[1]
    return(
        suppressMessages(corrr::correlate(df)) %>% 
            corrr::focus(., parentColumn) %>% 
            dplyr::filter(., parentColumn >= as.numeric(opt$corrThreshold)) %>% 
            dplyr::pull(., term)
    )
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

    # node$abundance <- hData[row,2:ncol(hData)]
    node$abundance <- hData[row,2:5]
    node$passedPrevalenceFilter <- rowSums(node$abundance != 0) > (NCOL(node$abundance)*filterPrevalence)
    node$passedMeanAbundanceFilter <- mean(unlist(node$abundance), trim = trim) > filterMeanAbundance
    node$winner <- FALSE    
    node$highlyCorrelated = FALSE
    node$id <- as.character(row)
}

hTree[["k__Archaea"]]$Do(function(node) node$competed <- {
    # handle no children
    if (length(node$children) == 0) {
        node$winner <- TRUE
        return(TRUE)
    }

    # build dataframe of parent and children
    # parent is row 1
    df <- node$abundance

    # handle descendant winners
    # take in depth limit
    descendantWinners <- getDescendantWinners(node, 5)
    if (length(descendantWinners) == 0) {
        node$winner <- TRUE
        return(TRUE)
    }
    for (i in 1:length(descendantWinners)) {
        df <- rbind(df, descendantWinners[[i]]$abundance)
    }

    # transpose
    transposed <- t(df)
    print(df)
    return()

    # standalone func
    # inputs: transposed dataframe, correlation threshold
    # mark correlated values
    # outputs: vector of correlated ids
    correlatedIDs <- calculateCorrelation(transposed, corrThreshold)

    # mark correlated in tree
    for (i in 1:length(winners)) {
        if (winners[[i]]$id %in% correlatedIDs) {
            winners[[i]]$winner <- FALSE
            winners[[i]]$highlyCorrelated <- TRUE
        }
    }
    # handle all children correlated, parent winner
    if (length(winners) == length(correlatedIDs)) {
        node$winner <- TRUE
        return(TRUE)
    }
    # drop from transposed data all correlated children

    # merge metadata??

    # standalone func
    # inputs: transposed dataframe, metadata, ...
    # random forest step
    # select factor
    # run model
    # outputs: vector of winner ids

    # mark winners/losers of parent and descendants
    
}, traversal = "post-order")

print(hTree[["k__Archaea"]], "level")
