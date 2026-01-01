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
    DIET_ML_TAGS="-t aoliver44/diet_ml:latest"

    # set platforms for max machine capability
fi

PLATFORM_FLAG="--platform linux/amd64,linux/arm64"

echo "Building taxaHFE containers for version $VERSION"

echo "Building individual command containers..."
# build args specifies which cmd/ script to use
docker buildx build --target taxa_hfe $PLATFORM_FLAG --build-arg version=$VERSION -t aoliver44/taxa_hfe:$VERSION $TAXA_HFE_TAGS --push .
docker buildx build --target taxa_hfe_ml $PLATFORM_FLAG --build-arg version=$VERSION -t aoliver44/taxa_hfe_ml:$VERSION $TAXA_HFE_ML_TAGS --push .
docker buildx build --target diet_ml $PLATFORM_FLAG --build-arg version=$VERSION -t aoliver44/diet_ml:$VERSION $DIET_ML_TAGS --push .

echo "Building the rstudio development image..."
docker buildx build --target rstudio $PLATFORM_FLAG --build-arg version=$VERSION -t aoliver44/taxa_hfe_rstudio:$VERSION $RSTUDIO_TAGS --push .
