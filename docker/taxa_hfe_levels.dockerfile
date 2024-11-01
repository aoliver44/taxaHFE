## base r image where we have already installed the packages
ARG BASE_IMAGE=aoliver44/taxa_hfe_base
FROM $BASE_IMAGE

# copy the desired command file
COPY cmd/*.R .

# set the data directory to the place where the data should be mounted
ENTRYPOINT [ "Rscript", "summarized_levels.R", "--data_dir", "/data"]