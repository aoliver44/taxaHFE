## main file for running taxaHFE
## Build: docker build --platform linux/amd64 -t aoliver44/taxa_hfe_base:{version} -t aoliver44/taxa_hfe_base:latest .
## once this is built, the individual containers for the different command can be built using this image
## Run: TODO

## base image to start with
FROM rocker/r-ver:4.5.2

## taxaHFE version, read in from `--build-arg version={}` in the docker build command
ARG version
ENV TAXA_HFE_VERSION=${version}

## RENV version
ENV RENV_VERSION=1.1.5

RUN apt-get update
RUN apt-get install -y libz-dev libxml2-dev libcurl4-openssl-dev libssl-dev libpng-dev python3

## install RENV the suggested way: https://rstudio.github.io/renv/articles/docker.html#creating-docker-images-with-renv
RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"
## install remotes bc its nice to have in the images esp for development
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"

WORKDIR /app

# copy in the renv lockfile and install all required packages
COPY renv.lock .
RUN R -e 'renv::consent(provided = TRUE)'
RUN R -e 'renv::restore()'

# copy in lib
COPY lib lib/
