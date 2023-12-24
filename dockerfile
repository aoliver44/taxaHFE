## Author: Andrew Oliver
## Version: aoliver44/taxaHFE:1.12
## Date: June 14, 2023

## base image to start with
FROM rocker/rstudio:4.2.0

## RENV version
ENV RENV_VERSION=0.16.0

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
COPY run_taxaHFEv2.R ./scripts/taxaHFE
COPY tree.R ./scripts/utilities/tree.R

ENV PATH="${PATH}:/scripts/"

# USER docker
