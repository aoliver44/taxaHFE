## develop on taxaHFE using rstudio
## Build: docker build --platform linux/amd64 -f rstudio.dockerfile -t aoliver44/taxa_hfe_rstudio:2.3 .
## Run: 
##  - docker run -p 8787:8787 -e PASSWORD={PASSWORD_HERE} -v `pwd`:/home/rstudio aoliver44/taxa_hfe_rstudio
##  - login at http://localhost:8787, user: rstudio, password: whatever was set in {PASSWORD_HERE}

## base r image where we have already installed the packages
ARG BASE_IMAGE=aoliver44/taxa_hfe_base:latest
FROM $BASE_IMAGE

## base rstudio image
FROM rocker/rstudio:4.2.3

## taxaHFE version, read in from `--build-arg version={}` in the docker build command
ARG version
ENV TAXA_HFE_VERSION=${version}

## RENV version
ENV RENV_VERSION=0.16.0

RUN apt-get update
RUN apt-get install -y libz-dev libxml2-dev python3

# grab packages installed in the base taxaHFE image
COPY --from=0 /usr/local/lib/R/site-library /usr/local/lib/R/site-library
COPY --from=0 /usr/local/lib/R/library /usr/local/lib/R/library
