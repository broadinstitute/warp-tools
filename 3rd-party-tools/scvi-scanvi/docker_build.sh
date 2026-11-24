#!/bin/bash
# Builds scvi-scanvi locally for testing. CI builds and pushes the image that WARP pins.
set -e
DOCKER_IMAGE_VERSION=1.7.3
SCVI_TOOLS_VERSION="1.5.1"  # the Dockerfile's ARG SCVI_TOOLS_VERSION default
TIMESTAMP=$(date +"%s")
IMAGE_TAG="$DOCKER_IMAGE_VERSION-$SCVI_TOOLS_VERSION-$TIMESTAMP"
docker build -t "scvi-scanvi:$IMAGE_TAG" "$(dirname "$0")"
