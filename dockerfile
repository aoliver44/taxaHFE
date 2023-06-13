## Author: Andrew Oliver
## Version: aoliver44/taxaHFE:1.10
## Date: May 3, 2023

## base image to start with
FROM rocker/r-base:4.2.0

## RENV version
ENV RENV_VERSION=0.17.3

RUN apt update
# install some things that R needs
RUN apt install -y libz-dev libxml2-dev

# install RENV, which will then install all R project packages
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

# should be in the same directory as this file
COPY renv.lock ./
RUN R -e 'renv::consent(provided = TRUE)'
RUN R -e 'renv::restore()'

# copy in scripts so they are part of container
COPY taxaHFE.R ./scripts/taxaHFE
COPY taxaHFE_functions.R ./scripts/utilities/taxaHFE_functions.R

ENV PATH="${PATH}:/scripts/"

USER docker
