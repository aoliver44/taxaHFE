set -e

# load version from env, falling back to file
# exit if netiher populate
if [ -z "${VERSION}" ]; 
    then VERSION=$(head -n 1 VERSION 2>/dev/null)
fi

if [ -z "${VERSION}" ]; 
    then echo "exiting, VERSION not set by env or file" && exit 1 
fi

echo "Building taxaHFE containers for version $VERSION"

echo "Building the base image and the rstudio image..."
docker build --platform linux/amd64 -f docker/base.dockerfile --build-arg version=$VERSION -t aoliver44/taxahfe_base:$VERSION -t aoliver44/taxahfe_base:latest .
docker build --platform linux/amd64 -f docker/rstudio.dockerfile --build-arg version=$VERSION -t aoliver44/taxahfe_rstudio:$VERSION -t aoliver44/taxahfe_rstudio:latest .

echo "Building individual command containers..."
# build args specifies which cmd/ script to use
docker build --platform linux/amd64 -f docker/taxahfe.dockerfile -t aoliver44/taxahfe:$VERSION -t aoliver44/taxahfe:latest .
docker build --platform linux/amd64 -f docker/taxahfe_ml.dockerfile -t aoliver44/taxahfe_ml:$VERSION -t aoliver44/taxahfe_ml:latest .
docker build --platform linux/amd64 -f docker/taxahfe_levels.dockerfile -t aoliver44/taxahfe_levels:$VERSION -t aoliver44/taxahfe_levels:latest .

# TODO docker login/push