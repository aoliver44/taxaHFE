FROM aoliver44/taxahfe_base:latest

# copy the desired command file
COPY cmd/*.R .

ENTRYPOINT [ "Rscript", "summarized_levels.R"]