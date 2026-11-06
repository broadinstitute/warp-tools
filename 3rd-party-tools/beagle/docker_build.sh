#!/bin/bash
set -e

# NOTE: This build script probably doesn't do what you hope it does.
# While it builds the image, it doesn't actually push it to GCR because
# the user running the script will most likely not have permissions to
# push to the GCR repo. The only utility of this script is to create a
# tag with the correct format and build the image locally. If you want to
# push to GCR, you'll have to do it using the GitHub action (Beagle CI)
# and specifying the correct tag (which you can get from running this script,
# or constructing manually).

# Update version when changes to Dockerfile are made
DOCKER_IMAGE_VERSION=1.1.0
TIMESTAMP=$(date +"%s")
DIR=$(cd $(dirname $0) && pwd)

# Registries and tags
GCR_URL="us.gcr.io/broad-gotc-prod/imputation-beagle"

# Beagle version
BEAGLE_VERSION="17Dec24.224"

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

    echo "building and pushing GCR Image - $GCR_URL:$IMAGE_TAG"

    # consider: add `--squash` when ready to productionize. https://docs.docker.com/reference/cli/docker/image/build/#squash
    docker build -t "$GCR_URL:$IMAGE_TAG" \
        --build-arg BEAGLE_VERSION="$BEAGLE_VERSION" \
        --no-cache $DIR
    docker push "$GCR_URL:$IMAGE_TAG"

    echo -e "$GCR_URL:$IMAGE_TAG" >> "$DIR/docker_versions.tsv"
    echo "done"
}

main "$@"
