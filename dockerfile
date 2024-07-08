## main file for running taxaHFE
## Build: docker build --platform linux/amd64 -t aoliver44/taxa_hfe:2.3 .
## Run: TODO

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

WORKDIR /app

# copy in the renv lockfile and install all required packages
COPY renv.lock .
RUN R -e 'renv::consent(provided = TRUE)'
RUN R -e 'renv::restore()'

# copy in scripts so they are part of container
COPY lib lib/
COPY models models/
COPY *.R .

ENTRYPOINT ["Rscript", "leakfree_taxaHFE.R"]
