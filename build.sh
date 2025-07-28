set -e

if [ $1 == '--dev' ]; then
    VERSION="dev"
else
    # if not a dev build, load VERSION and also tag latest
    # load version from env, falling back to file
    # exit if netiher populate
    if [ -z "${VERSION}" ]; then 
        VERSION=$(head -n 1 VERSION 2>/dev/null)
    fi

    if [ -z "${VERSION}" ]; then
        echo "exiting, VERSION not set by env or file" && exit 1 
    fi

    BASE_TAGS="-t aoliver44/taxa_hfe_base:latest"
    RSTUDIO_TAGS="-t aoliver44/taxa_hfe_rstudio:latest"
    TAXA_HFE_TAGS="-t aoliver44/taxa_hfe:latest"
    TAXA_HFE_ML_TAGS="-t aoliver44/taxa_hfe_ml:latest"
fi

BUILD_ARGS="--build-arg version=$VERSION --build-arg BASE_IMAGE=aoliver44/taxa_hfe_base:$VERSION"

echo "Building taxaHFE containers for version $VERSION"

echo "Building the base image and the rstudio image..."
docker build --platform linux/amd64 -f docker/base.dockerfile $BUILD_ARGS -t aoliver44/taxa_hfe_base:$VERSION $BASE_TAGS .
docker build --platform linux/amd64 -f docker/rstudio.dockerfile $BUILD_ARGS -t aoliver44/taxa_hfe_rstudio:$VERSION $RSTUDIO_TAGS .

echo "Building individual command containers..."
# build args specifies which cmd/ script to use
docker build --platform linux/amd64 -f docker/taxa_hfe.dockerfile $BUILD_ARGS -t aoliver44/taxa_hfe:$VERSION $TAXA_HFE_TAGS .
docker build --platform linux/amd64 -f docker/taxa_hfe_ml.dockerfile $BUILD_ARGS -t aoliver44/taxa_hfe_ml:$VERSION $TAXA_HFE_ML_TAGS .

# TODO docker login/push