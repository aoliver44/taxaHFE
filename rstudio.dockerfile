## main file for running taxaHFE
## Build: docker build --platform linux/amd64 -f rstudio.dockerfile -t aoliver44/taxa_hfe_rstudio:2.3 .
## Run: 
##  - docker run -p 8787:8787 -e PASSWORD={PASSWORD_HERE} -v `pwd`:/home/rstudio test_taxa_hfe_rstudio
##  - login at http://localhost:8787, user: rstudio, password: whatever was set in {PASSWORD_HERE}

## base image to start with
FROM rocker/rstudio:4.2.3

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

# copy in the renv lockfile and install all required packages
COPY renv.lock .
RUN R -e 'renv::consent(provided = TRUE)'
RUN R -e 'renv::restore()'
