#remotes::install_github(repo = "cran/htree")
library(htree)

metadata <- read_in_metadata(input = "example_inputs/metadata_time.txt", subject_identifier = "subject_id", label = "Intervention")
hData <- read_in_hierarchical_data(input = "example_inputs/microbiome_time.txt", metadata = metadata, cores = 4)

metadata <- metadata %>% 
  group_by(subject) %>% 
  mutate(time_order = order(subject, week_num))

hTree <- build_tree(hData,
                    filter_prevalence = 0.1,
                    filter_mean_abundance = 0.01
)



## setting classify=TRUE builds classification tree (gini-impurity node splitting)
practice_data <- read.csv("example_inputs/taxaHFE_time_input.csv")
practice_data <- practice_data %>% drop_na() %>% select(., -subject_id)
#historical_predictors=match(c("stress","illness"),names(mscm))
concurrent_predictors=c(3,4,6,7:12)
historical_predictors=c(7:12)
control=list(se=TRUE, classify=TRUE,vc=concurrent_predictors,vh=historical_predictors)
ff=hrf(x=practice_data,id=practice_data$subject,time=practice_data$week_num,yindx="feature_of_interest")

ff=hrf(x=merged_data,id=merged_data$individual,time=merged_data$time,yindx="feature_of_interest", control=list(se=TRUE, classify=TRUE, B = 2, R = 3))

vi=varimp_hrf(ff, nperm = 40, parallel = FALSE)
vi %>% dplyr::select(1, 4) %>% dplyr::rename(., "taxa" = 1, "importance" = 2)

ranger::ranger(as.factor(feature_of_interest) ~ ., data = practice_data, importance = "impurity_corrected", sample.fraction = 1, replace = TRUE, num.threads = 4)$variable.importance %>%
  as.data.frame() %>%
  dplyr::rename(., "importance" = ".") %>%
  tibble::rownames_to_column(var = "taxa")


data(cd4)
data("mscm")
control=list(se=TRUE)
ff=hrf(x=cd4,id=cd4$id,time=cd4$time,yindx="count",control=control)

vi=varimp_hrf(ff)
vi
# -- partial dependence for top 4 predictors (with +/-2 SE estimates) 
par(mfrow=c(2,2))
for(k in 1:4)
  pd=partdep_hrf(ff,xindx=as.character(vi$Predictor[k]))
par(mfrow=c(1,1))

plot(1:length(ff$error),ff$error,xlab="forest size",ylab="oob mse",type="l")

## by default, the number of delta values (parameter 'eta_1' above) is 20
## can set this using 'ndelta'  
control$ndelta=50

control$se=FALSE # -- turning off bootstrapping .. 
ff=hrf(x=cd4,id=cd4$id,time=cd4$time,yindx="count",control=control)
points(1:length(ff$error),ff$error,type="l",lty=2)

# the grid of delta values 
ff$control$delta

historical_predictors=match(c("stress","illness"),names(mscm))
concurrent_predictors=which(names(mscm)!="stress")


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
    control=list(se=TRUE, classify=TRUE)
  } else {
    response_formula <- as.formula(paste("as.numeric(", feature_of_interest, ") ~ .", sep = ""))
    control=list(se=TRUE, classify=FALSE)
  }
  
  # progress bar for the final rf competition
  # will only be incremented/shown if parent_descendent_competition == FALSE
  pb <- progress::progress_bar$new(format = " Competing final winners [:bar] :percent in :elapsed", total = nperm, clear = FALSE, width = 60)
  
  # run ranger, setting parameters such as
  # random seed
  # num.threads number of threads to five ranger
  run_rf <- function(seed) {
    if (!parent_descendent_competition) pb$tick()
    
    if (random_effects) {
      
      re_rf <- hrf(x=merged_data,id=merged_data$individual,time=merged_data$day,yindx="feature_of_interest", control = control)
      re_vi=varimp_hrf(re_rf, nperm = nperm, parallel = ncores)
      re_vi %>% dplyr::select(1, 4) %>% dplyr::rename(., "taxa" = 1, "importance" = 2)
      
    } else {
      ranger::ranger(response_formula, data = merged_data, importance = "impurity_corrected", seed = seed, sample.fraction = 1, replace = TRUE, num.threads = ncores)$variable.importance %>%
        as.data.frame() %>%
        dplyr::rename(., "importance" = ".") %>%
        tibble::rownames_to_column(var = "taxa")
    }
  }
  
  # run the above function across nperm random seeds and average the vip scores
  model_importance <- purrr::map_df(sample(1:1000000, nperm), run_rf) %>%
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

library(dplyr)
count_reads <- read.delim("example_inputs/count_reads.txt", sep = "\t")
genome_taxonomy <- read.delim("example_inputs/genome_taxonomy.txt", sep = "\t")
genome_info <- read.delim("example_inputs/genome_info.txt", sep = "\t")
metadata <- read.delim("example_inputs/fiber_metadata.txt", sep = "\t")

genome_info <- genome_info %>% select(., genome_id, species_id) %>% filter(., !duplicated(species_id))
genome_taxonomy <- genome_taxonomy %>% select(., genome_id, kingdom:species)
genome_taxonomy$genome_id <- gsub(pattern = "00$", replacement = "", x = genome_taxonomy$genome_id)

count_reads <- genome_info %>% left_join(., count_reads)
count_reads <- merge(genome_taxonomy, count_reads, by = "genome_id")

count_reads$clade_name <- paste0("k__", count_reads$kingdom, "|p__", count_reads$phylum, "|c__", count_reads$class, "|o__", count_reads$order, "|f__", count_reads$family, "|g__", count_reads$genus, "|s__", count_reads$species_id)
count_reads <- count_reads %>% 
  relocate(., clade_name) %>% 
  select(., -c(2:10)) %>% 
  tibble::column_to_rownames(., var = "clade_name")

count_reads$rowsums <- rowSums(count_reads)
count_reads <- count_reads %>% filter(., rowsums > 0)
library(EcolUtils)
library(vegan)
set.seed(seed = 999)
#barplot(sort(rowSums(midas)))
#abline(h = mean(rowSums(midas)), col = "Red")
#abline(h = median(rowSums(midas)), col = "Blue")
count_reads_rare <- rrarefy.perm(t(count_reads), sample = 900, n = 100, round.out = T)
count_reads_rare <- count_reads_rare[rowSums(count_reads_rare) >= 900-(900*.1), 
                               colSums(count_reads_rare) >= 1]

count_reads_rare <- count_reads_rare %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., var = "clade_name")

write.table(x = count_reads_rare, file = "example_inputs/microbiome_time.txt", quote = F, row.names = F, sep = "\t")

microbiome <- read.delim("example_inputs/microbiome_time.txt", sep = "\t")
metadata <- read.delim("example_inputs/metadata_time.txt", sep = "\t")

## need to fix "alistipes o." by hand or make some code to do it