# TODO - is this needed?

#!/bin/bash
set -e

# NOTE: This script builds the image locally with the correct tag format
# but doesn't automatically push to GCR. To push to GCR, use the GitHub action
# workflow and specifying the correct tag (which you can get from running
# this script, or constructing manually).

# Update version when changes to Dockerfile are made
DOCKER_IMAGE_VERSION=0.0.1
TIMESTAMP=$(date +"%s")
DIR=$(cd $(dirname $0) && pwd)

# Registries and tags
GCR_URL="us.gcr.io/broad-gotc-prod/glimpse"

# GLIMPSE commit hash - update this when you want to build from a different commit
GLIMPSE_COMMIT_HASH="cf5822c326d9e4b0e4d2a273fd3497aa3d62f7f0"

# Necessary tools and help text
TOOLS=(docker gcloud)
HELP="$(basename "$0") [-h|--help] [-c|--commit] [-t|--tools] -- script to build the GLIMPSE image and push to GCR

where:
    -h|--help Show help text
    -c|--commit Commit hash of GLIMPSE to use (default: GLIMPSE_COMMIT_HASH=${GLIMPSE_COMMIT_HASH})
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
        -c|--commit)
        GLIMPSE_COMMIT_HASH="$2"
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

    # Use short commit hash (first 7 characters) for the tag
    SHORT_COMMIT_HASH="${GLIMPSE_COMMIT_HASH:0:7}"
    IMAGE_TAG="$DOCKER_IMAGE_VERSION-$SHORT_COMMIT_HASH-$TIMESTAMP"

    echo "building and pushing GCR Image - $GCR_URL:$IMAGE_TAG"

    # consider: add `--squash` when ready to productionize. https://docs.docker.com/reference/cli/docker/image/build/#squash
    docker build -t "$GCR_URL:$IMAGE_TAG" \
        --build-arg GLIMPSE_COMMIT_HASH="$GLIMPSE_COMMIT_HASH" \
        --no-cache $DIR
    docker push "$GCR_URL:$IMAGE_TAG"

    echo -e "$GCR_URL:$IMAGE_TAG" >> "$DIR/docker_versions.tsv"
    echo "done"
}

main "$@"
