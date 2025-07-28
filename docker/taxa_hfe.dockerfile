## base r image where we have already installed the packages
ARG BASE_IMAGE=aoliver44/taxa_hfe_base
FROM $BASE_IMAGE

# copy the desired command file
COPY cmd/*.R .

# set the data directory to the place where the data will be mounted
ENTRYPOINT [ "Rscript", "taxa_hfe.R", "--data_dir", "/data"]