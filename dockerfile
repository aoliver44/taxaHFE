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

# install RENV, which will then install all R project packages
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

# should be in the same directocat ry as this file
COPY renv.lock ./
RUN R -e 'renv::consent(provided = TRUE)'
RUN R -e 'renv::restore()'

# copy in scripts so they are part of container
COPY ./taxaHFE-SHAP/leakfree_taxaHFE.R ./scripts/taxaHFE-ML
COPY ./taxaHFE-SHAP/dietML.R ./scripts/dietML.R
COPY ./tree.R ./scripts/tree.R
COPY ./run_taxaHFEv2.R ./scripts/taxaHFE
COPY ./taxaHFE-SHAP/models/dietML_ranger_tidy.R ./scripts/models/dietML_ranger_tidy.R
COPY ./taxaHFE-SHAP/models/dietML_null_tidy.R ./scripts/models/dietML_null_tidy.R
COPY ./taxaHFE-SHAP/utilities/shap_figures.R ./scripts/utilities/shap_figures.R
COPY ./taxaHFE-SHAP/utilities/vip_basic.R ./scripts/utilities/vip_basic.R

ENV PATH="${PATH}:/scripts/"


