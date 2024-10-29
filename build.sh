set -e

# load version from env, falling back to file
# exit if netiher populate
if [ -z "${VERSION}" ]; then 
    VERSION=$(head -n 1 VERSION 2>/dev/null)
fi

if [ -z "${VERSION}" ]; then
    echo "exiting, VERSION not set by env or file" && exit 1 
fi

BASE_TAGS="--build-arg version=$VERSION -t aoliver44/taxahfe_base:$VERSION -t aoliver44/taxahfe_base:latest"
RSTUDIO_TAGS="--build-arg version=$VERSION -t aoliver44/taxahfe_rstudio:$VERSION -t aoliver44/taxahfe_rstudio:latest"
TAXAHFE_TAGS="-t aoliver44/taxahfe:$VERSION -t aoliver44/taxahfe:latest"
TAXAHFE_ML_TAGS="-t aoliver44/taxahfe_ml:$VERSION -t aoliver44/taxahfe_ml:latest"
TAXAHFE_LEVELS_TAGS="-t aoliver44/taxahfe_levels:$VERSION -t aoliver44/taxahfe_levels:latest"

if [ $1 == '--dev' ]; then
    VERSION="dev"
    BASE_TAGS="--build-arg version=$VERSION -t aoliver44/taxahfe_base:$VERSION"
    RSTUDIO_TAGS="--build-arg version=$VERSION -t aoliver44/taxahfe_rstudio:$VERSION"
    TAXAHFE_TAGS="-t aoliver44/taxahfe:$VERSION"
    TAXAHFE_ML_TAGS="-t aoliver44/taxahfe_ml:$VERSION"
    TAXAHFE_LEVELS_TAGS="-t aoliver44/taxahfe_levels:$VERSION"
fi

echo "Building taxaHFE containers for version $VERSION"

echo "Building the base image and the rstudio image..."
docker build --platform linux/amd64 -f docker/base.dockerfile $BASE_ARGS .
docker build --platform linux/amd64 -f docker/rstudio.dockerfile $RSTUDIO_TAGS .

echo "Building individual command containers..."
# build args specifies which cmd/ script to use
docker build --platform linux/amd64 -f docker/taxahfe.dockerfile $TAXAHFE_TAGS .
docker build --platform linux/amd64 -f docker/taxahfe_ml.dockerfile $TAXAHFE_ML_TAGS .
docker build --platform linux/amd64 -f docker/taxahfe_levels.dockerfile $TAXAHFE_LEVELS_TAGS .

# TODO docker login/push