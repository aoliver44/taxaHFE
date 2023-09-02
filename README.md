
 # **taxaHFE**
 
 A program to perform hierarchical feature engineering on data with taxonomic organization (i.e., microbiome data, dietary data)

## Table of Contents
- [Download taxaHFE](https://github.com/aoliver44/taxaHFE#download-taxahfe)
- [Quickstart](https://github.com/aoliver44/taxaHFE#quickstart)
- [Flag information](https://github.com/aoliver44/taxaHFE#information-about-the-flags)
- [About](https://github.com/aoliver44/taxaHFE#about)
- [Contribute](https://github.com/aoliver44/taxaHFE#contribute)



-----------------------------



## **Outline of taxaHFE**


![Outline of taxaHFE algorithm](algorithm_outline_github.png "Outline of taxaHFE algorithm")

</br>

------------------------------
## Download taxaHFE


Option 1: The easiest way to get started is pulling the docker image. Please [install docker](https://www.docker.com/) you go this route. 

```
docker pull aoliver44/taxa_hfe:2.0
```

Option 2: Alternatively, you can pull this image using Singularity:

```
singularity pull taxaHFE.sif docker://aoliver44/taxa_hfe:2.0
```

Option 3: Finally, it's possible to build the image yourself:

1. Download the dockerfile and renv.lock file from github
2. Navigate to the directory with these files
3. Run the command:

```
docker build -t taxa_hfe:2.0 .
```
</br>

------------------------------
## Quickstart

Option 1: Run taxaHFE with **YOUR** data:
1. Navigate to the directory containing your data, and start the docker image!
```
docker run --rm -it -v `pwd`:/home/docker -w /home/docker aoliver44/taxa_hfe:2.0 bash

## or with singularity
singularity run -W `pwd` --bind `pwd`:/home/docker taxaHFE.sif bash
```
2. Run taxaHFE
```
taxaHFE.R [options] <METADATA> <DATA> <OUTPUT>
```
OR

Option 2: Run taxaHFE on **EXAMPLE** data provided:

```
## STEP 1: CLONE THE REPOSITORY
git clone https://github.com/aoliver44/taxaHFE.git && cd taxaHFE

## STEP 2: RUN THE CONTAINER
docker run --rm -it -v `pwd`:/home/docker -w /home/docker aoliver44/taxa_hfe:2.0 bash

## STEP 3: RUN TAXAHFE 
taxaHFE --subject_identifier Sample --label Category --lowest_level 3 --ncores 2 --seed 42 /home/docker/example_inputs/metadata.txt /home/docker/example_inputs/microbiome_data.txt /home/docker/example_inputs/output.csv
```

</br>

------------------------------
## Information about the flags


```
Hierarchical feature engineering (HFE) for the reduction of features with respects to a factor or regressor

usage:    
taxaHFE.R [options] <METADATA> <DATA> <OUTPUT>    

Options:   
 -h --help  Show this screen.    
 -v --version  Show version.    
 -s --subject_identifier metadata column name containing subject IDs [default: subject_id]    
 -l --label metadata column name of interest for ML [default: cluster]    
 -t --feature_type is the ML label a factor or numeric [default: factor]    
 -f --sample_fraction only let rf see a fraction of total data [default: 1]    
 -a --abundance feature abundance filter [default: 0.0001]    
 -p --prevalence feature prevalence filter [default: 0.01]    
 -L --lowest_level most general level allowed to compete [default: 2]    
 -m --max_depth how many hierarchical levels should be allowed to compete [default: 1000]    
 -c --cor_level initial pearson correlation filter [default: 0.95]    
 -w --write_old_files write individual level files and old HFE files [default: TRUE]    
 -n --ncores number of cpu cores to use [default: 2]   
 --seed set a random seed, default is to use system time 
 
 Arguments:    
 
 METADATA path to metadata input (txt | tsv | csv)
 DATA path to input file from hierarchical data (i.e. hData data) (txt | tsv | csv)
 OUTPUT output file name (csv)
```

--subject_identifier: this is a column that identifies the sample or subject ID in the input metadata. All subjectIDs should be unique. They will be coerced to unique values (and simplified snake_case alpha-numerics) using ```janitor::make_clean_names()```

--label: the name of the column in your input metadata that you are trying to predict with HFE. Can be a factor or continous. 

--feature_type: is the label a factor or a continous variable (options: factor or numeric)?

--subsample: a decimal value for performing stratified subsampling of factor type data. This behavior is to help protect against data-leakage.

--abundance: a per-feature abundance filter. This filter calculates an outlier-resistant mean (trimming the top and bottom 2% of data) of the feature's abundance. If the average abundance across the middle 96% of samples is above this value, the feature is kept. Note, if your sampling effort is not standardized in some way (e.g. relative abundance), this filter may produce undesirable behavior. To turn this filter off, set its value to 0 (or the minimum value in your dataset). The default behavior is to filter out features below a mean abundance of 0.0001; however, this assumes the feature abundances exist on a scale from 0-1. 

--prevalence: a per-feature prevalence filter. This filter sets the number of non-zero occurances desired for features. The default behavior is if the feature is 99% zeros, it will be dropped from further analysis. This filter is also somewhat sensitive to sampling depth, as samples with greater sampling depth will likely find rarer features.

--lowest_level: The lowest level for which to compete in a taxaHFE competition. To better understand this parameter, consider a microbiome competition as an example: Each feature contains some version of taxonomic levels from general -> specific (kingdom, phylum, class, order, family, genus, species). taxaHFE adds one additional level "below" kingdom, called taxa_tree (a somewhat meaningless root repersenting the sum-total abundance per sample). Setting ```--lowest_level 1``` allows taxaHFE to take competitions all the way to taxa_tree, potentially allowing it to be the only feature selected (if, for instance, the differences in sum-total abundances are the most informative feature with respects to your metadata label). The default behavior is to set ```--lowest_level 2```, which would stop the competitions at the kingdom level in this example. Sometimes, this behavior will result in a similar result desicribed above (i.e., the kingdom Bacteria is selected as the sole winner). IF your interest is what features are informative within your most general grouping (i.e., which bacteria, archaea, etc.), then considering setting this value to ```--lowest_level 3```.

--max_depth: how deep should a child be allowed to compete? In version 1 of this program, max_depth was effectively 1, which meant that a child was only allowed to compete against their parent, but NOT their grandparent. The default behavior in the current version is ```--max_depth 1000```, which means a child is allowed to compete against their parent and 1000 hierarchical levels beyond their parent (i.e., great-great-great...grandparent). If they are an informative feature, they are allowed to keep competing.

--cor_level: what initial correlation threshold (Pearson) to use when comparing child to parent. We use a high threshold (0.95) and encourage this threashold to stay high. We are really after *redundant* features with this step, were are not trying to institute a deep correlation filter.

--write_old_files: should files summarized at each taxa level be written to file? The old HFE program files are written to for use in the Oudah et al. algorithm

--ncores: number of cores to let the random forest use. Not a huge spead up, but if you have the cores available, it can't hurt.

--seed: the default behavior is to use ```Sys.time()``` to generate a random seed each time taxaHFE is run. If you set it to a number, it will likely return the same results across repeated runs (though this assumption has not been thoroughly tested).

[METADATA]: A **full path** the file that contains the metadata column you wish to predict with your hierarchical data. This file should contain BOTH your subject_identifier and your metadata label

[DATA]: A **full path** to your taxonomic or hierarchical feature set. Columns should be your subject_identifier, plus one column labeled clade_name. Each feature should have the levels seperated using a pipe seperator (```"|"```).

[OUTPUT]: A **full path** to a file which will be the main output file. The other output files will parse this path. 


**OUTPUTS**

**output_level_[1,2,3...].csv:** summarized files at each taxonomic level (if write_old_files = TRUE)

**output.txt:** taxaHFE processed dataset, with super filter

**output_nosf.txt:** taxaHFE processed dataset, without super filter

**output_old_hfe_label.txt:** Label data for use in the Oudah algorithm (if write_old_files = TRUE)

**output_old_hfe_otu.txt:** OTU data for use in the Oudah algorithm (if write_old_files = TRUE)

**output_old_hfe_taxa.txt:** Taxa data for use in the Oudah algorithm (if write_old_files = TRUE)

</br>

------------------------------
## About

We developed software, called taxaHFE (Hierarchical Feature Engineering), which works by first considering the pairwise correlation structure between a taxon and its descendants to prune descendants above a correlation threshold. Next it permutes a random forest on the taxon and remaining descendants to determine how important each is at explaining an intervention or clinical covariate. If, on average, the taxon is the most important feature in the model, the descendants are dropped, otherwise only the descendants more important than the taxon are kept. Last, an optional final filter step considers all features remaining, and again permutes a random forest. Any features which are either below the average importance of all remaining features or have a negative or zero average importance are dropped.  

</br>

------------------------------
## Contribute

Feel free to raise an issue, contribute with a pull request, or reach out!





