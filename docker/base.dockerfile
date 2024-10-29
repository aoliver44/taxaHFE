## main file for running taxaHFE
## Build: docker build --platform linux/amd64 -t aoliver44/taxahfe_base:{version} -t aoliver44/taxahfe_base:latest .
## once this is built, the individual containers for the different command can be built using this image
## Run: TODO

## base image to start with
FROM rocker/r-ver:4.2.3

## taxaHFE version, read in from `--build-arg version={}` in the docker build command
ARG version
ENV TAXAHFE_VERSION=${version}

## RENV version
ENV RENV_VERSION=0.16.0

RUN apt-get update
RUN apt-get install -y libz-dev libxml2-dev python3

RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

WORKDIR /app

# copy in the renv lockfile and install all required packages
COPY renv.lock .
RUN R -e 'renv::consent(provided = TRUE)'
RUN R -e 'renv::restore()'

# copy in lib
COPY lib lib/
