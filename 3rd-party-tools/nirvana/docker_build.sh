#!/bin/bash
set -e

# Update version when changes to Dockerfile are made
DOCKER_IMAGE_VERSION="v1.0.0"
TIMESTAMP=$(date +"%s")
DIR=$(cd "$(dirname "$0")" && pwd)

# Registry and tags
GCR_URL="us.gcr.io/broad-gotc-prod/nirvana"

# Nirvana Version and Configuration
NIRVANA_VERSION="v3.18.1"
DOTNET_CONFIGURATION="Release"

VERSION_NUMBER=$(echo "$NIRVANA_VERSION" | grep -Eo '[0-9]+([.][0-9]+)')
IMAGE_TAG="$DOCKER_IMAGE_VERSION-$VERSION_NUMBER-$TIMESTAMP"

# Necessary tools and help text
TOOLS=("docker")
HELP="$(basename "$0") [-h|--help] [-t|tools] -- script to build the Nirvana image and push to registry

where:
    -h|--help Show help text
    -t|--tools Show tools needed to run script
    "

function main(){
    for t in "${TOOLS[@]}"; do which "$t" >/dev/null || ok=no; done
    if [[ $ok == no ]]; then
        echo "Missing one of the following tools: "
        for t in "${TOOLS[@]}"; do echo "$t"; done
        exit 1
    fi

    while [[ $# -gt 0 ]]
    do 
        key="$1"
        case $key in
            -h|--help)
                echo "$HELP"
                exit 0
                ;;
            -t|--tools)
                for t in "${TOOLS[@]}"; do echo "$t"; done
                exit 0
                ;;
            *)
                shift
                ;;
        esac
    done

    # Build Nirvana Docker image
    echo "Building Nirvana Docker image... $GCR_URL:$IMAGE_TAG"
    docker build --no-cache -t "$GCR_URL:$IMAGE_TAG" \
        --build-arg NIRVANA_VERSION="$NIRVANA_VERSION" \
        --build-arg DOTNET_CONFIGURATION="$DOTNET_CONFIGURATION" \
        "$DIR"

    # Push to registry
    echo "Pushing Nirvana Docker image to registry..."
    docker push "$GCR_URL:$IMAGE_TAG"

    echo "$GCR_URL:$IMAGE_TAG" >> "$DIR/docker_versions.tsv"
    echo "Done"
}

main "$@"
