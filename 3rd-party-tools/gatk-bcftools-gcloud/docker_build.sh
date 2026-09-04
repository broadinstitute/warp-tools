#!/bin/bash
set -e

# Update version when changes to Dockerfile are made
DOCKER_IMAGE_VERSION=1.0.0
TIMESTAMP=$(date +"%s")
DIR=$(cd $(dirname $0) && pwd)

# Registries and tags
GCR_URL="us.gcr.io/broad-gotc-prod/gatk-bcftools-gcloud"

# GATK4 version
GATK4_VERSION="4.2.6.1"

# BCFTools version
BCFTOOLS_VERSION="1.24"

# Necessary tools and help text
TOOLS=(docker gcloud)
HELP="$(basename "$0") [-h|--help] [-v4|--version4] [-b|--bcf] [-t|--tools] -- script to build the GATK/BCFtools/gcloud image and push to GCR

where:
    -h|--help Show help text
    -v4|--version4 Version of GATK4 to use (default: GATK4=$GATK4_VERSION)
    -b|--bcf Version of BCFTOOLS to use (default: BCFTOOLS_VERSION=${BCFTOOLS_VERSION})
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
        -v4|--version4)
        GATK4_VERSION="$2"
        shift
        shift
        ;;
        -b|--bcf)
        BCFTOOLS_VERSION="$2"
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

    IMAGE_TAG="$DOCKER_IMAGE_VERSION-$GATK4_VERSION-$BCFTOOLS_VERSION-$TIMESTAMP"

    echo "building and pushing GCR Image - $GCR_URL:$IMAGE_TAG"
    docker build -t "$GCR_URL:$IMAGE_TAG" \
        --build-arg GATK4_VERSION="$GATK4_VERSION" \
        --build-arg BCFTOOLS_VERSION="$BCFTOOLS_VERSION" \
        --no-cache $DIR
    docker push "$GCR_URL:$IMAGE_TAG"

    echo -e "$GCR_URL:$IMAGE_TAG" >> "$DIR/docker_versions.tsv"
    echo "done"
}

main "$@"
