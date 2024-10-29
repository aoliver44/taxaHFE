## base r image where we have already installed the packages
ARG BASE_IMAGE=aoliver44/taxa_hfe_base
FROM $BASE_IMAGE

# copy the desired command file
COPY cmd/*.R .

ENTRYPOINT [ "Rscript", "summarized_levels.R"]