## Author: Andrew Oliver
## Version: aoliver44/leakage_free_taxaHFE:latest
## Date: Jan 25, 2023

## base image to start with
FROM rocker/r-ver:4.2.3

## RENV version
ENV RENV_VERSION=0.16.0

RUN apt update
# install some things that R needs
RUN apt install -y libz-dev libxml2-dev
RUN apt-get update
RUN apt-get install -y python3

# install RENV, which will then install all R project packages
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

# should be in the same directocat ry as this file
COPY renv.lock ./
RUN R -e 'renv::consent(provided = TRUE)'
RUN R -e 'renv::restore()'

# copy in scripts so they are part of container
COPY ./taxaHFE-ML/leakfree_taxaHFE.R ./scripts/taxaHFE-ML
COPY ./taxaHFE-ML/dietML.R ./scripts/dietML.R
COPY ./tree.R ./scripts/tree.R
COPY ./options.R ./scripts/options.R
COPY ./run_taxaHFEv2.R ./scripts/taxaHFE
COPY ./taxaHFE-ML/models/dietML_ranger_tidy.R ./scripts/models/dietML_ranger_tidy.R
COPY ./taxaHFE-ML/models/dietML_null_tidy.R ./scripts/models/dietML_null_tidy.R
COPY ./taxaHFE-ML/utilities/shap_figures.R ./scripts/utilities/shap_figures.R
COPY ./taxaHFE-ML/utilities/vip_basic.R ./scripts/utilities/vip_basic.R

ENV PATH="${PATH}:/scripts/"
