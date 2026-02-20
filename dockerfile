## Dockerfile for taxahfe
# this file has multi-stage builds for the cmd/ files and the rstudio deelopment environment

################################################
# base image to start for all builds

FROM rocker/r-ver:4.5.2 AS base

# taxaHFE version, read in from `--build-arg version={}` in the docker build command
ARG version
ENV TAXA_HFE_VERSION=${version}

# RENV version
ENV RENV_VERSION=1.1.5

RUN apt-get update
RUN apt-get install -y libz-dev libxml2-dev libcurl4-openssl-dev libssl-dev libpng-dev python3
RUN apt-get install -y pigz

# install RENV the suggested way: https://rstudio.github.io/renv/articles/docker.html#creating-docker-images-with-renv
RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"

WORKDIR /app

# copy in the renv lockfile and install all required packages
COPY renv.lock .
RUN R -e 'renv::consent(provided = TRUE)'
RUN R -e 'renv::restore()'

# copy in lib
COPY lib lib/
COPY cmd/*.R .

################################################
# cmd directory builds
# each one uses the base image and adds the correct entrypoint 

# taxa_hfe
FROM base AS taxa_hfe
ENTRYPOINT [ "Rscript", "taxa_hfe.R", "--data_dir", "/data"]

# taxa_hfe_ml
FROM base AS taxa_hfe_ml
ENTRYPOINT [ "Rscript", "taxa_hfe_ml.R", "--data_dir", "/data"]

# diet_ml
FROM base AS diet_ml
ENTRYPOINT [ "Rscript", "diet_ml.R", "--data_dir", "/data"]

################################################
# rstudio build for development
# uses the rstudio image but takes the installed r packages from the base image

FROM rocker/rstudio:4.5.2 AS rstudio

# taxaHFE version, read in from `--build-arg version={}` in the docker build command
ARG version
ENV TAXA_HFE_VERSION=${version}

# RENV version
ENV RENV_VERSION=1.1.5

RUN apt-get update
RUN apt-get install -y libz-dev libxml2-dev python3
RUN apt-get install -y pigz

# grab packages installed in the base taxaHFE image
COPY --from=base /usr/local/lib/R/site-library /usr/local/lib/R/site-library
COPY --from=base /usr/local/lib/R/library /usr/local/lib/R/library

# install remotes bc its nice to have for development
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
