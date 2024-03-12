#!/bin/bash
set -e

# Update version when changes to Dockerfile are made
DOCKER_IMAGE_VERSION=0.0.1
TIMESTAMP="wip-temp-20240301"   # $(date +"%s")
DIR=$(cd $(dirname $0) && pwd)

# Registries and tags
# GCR_URL="us.gcr.io/broad-gotc-prod/imputation-beagle"

# GAR setup
GAR_REGION="us-central1"
GAR_PROJECT="morgan-fieldeng-gcp"
GAR_REPOSITORY="imputation-beagle-development"
GAR_IMAGE="imputation-beagle"
GAR_URL="${GAR_REGION}-docker.pkg.dev/${GAR_PROJECT}/${GAR_REPOSITORY}/${GAR_IMAGE}"

# Beagle version
BEAGLE_VERSION="01Mar24.d36"

# Necessary tools and help text
TOOLS=(docker gcloud)
HELP="$(basename "$0") [-h|--help] [-b|--beagle] [-t|--tools] -- script to build the Imputation Beagle image and push to GAR

where:
    -h|--help Show help text
    -b|--beagle Version of Beagle to use (default: BEAGLE_VERSION=${BEAGLE_VERSION})
    -t|--tools Show tools needed to run script
    "

function main(){
    for t in "${TOOLS[@]}"; do which $t >/dev/null || ok=no; done
        if [[ $ok == no ]]; then
            echo "Missing one of the following tools: "
            for t in "${TOOLS[@]}"; do echo "$t"; done
            exit 1
        fi

    while [[ $# -gt 0 ]]
    do 
    key="$1"
    case $key in
        -b|--beagle)
        BEAGLE_VERSION="$2"
        shift
        shift
        ;;
        -h|--help)
        echo "$HELP"
        exit 0
        ;;
        -t|--tools)
        for t in "${TOOLS[@]}"; do echo $t; done
        exit 0
        ;;
        *)
        shift
        ;;
    esac
    done

    IMAGE_TAG="$DOCKER_IMAGE_VERSION-$BEAGLE_VERSION-$TIMESTAMP"

    echo "building and pushing GCR Image - $GAR_URL:$IMAGE_TAG"
    docker build -t "$GAR_URL:$IMAGE_TAG" \
        --build-arg BEAGLE_VERSION="$BEAGLE_VERSION" \
        $DIR   
        # --no-cache $DIR\
    docker push "$GAR_URL:$IMAGE_TAG"

    echo -e "$GAR_URL:$IMAGE_TAG" >> "$DIR/docker_versions.tsv"
    echo "done"
}

main "$@"