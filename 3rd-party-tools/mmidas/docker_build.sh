#!/bin/bash
set -e

# Update version when changes to Dockerfile are made
DOCKER_IMAGE_VERSION=1.0.0
TIMESTAMP=$(date +"%s")
DIR=$(cd "$(dirname "$0")" && pwd)

# Registry
GCR_URL="us.gcr.io/broad-gotc-prod/mmidas"

# MMIDAS version (matches setup.py / pyproject.toml in the local source tree)
MMIDAS_VERSION="0.1.0"

# Caching: ON by default for faster iterative development.
# Pass --no-cache for a fully clean rebuild when needed.
CACHING="ON"

# Local MMIDAS source directory — resolved relative to this script's location.
# Assumes workspace layout:
#   <root>/MMIDAS/             ← source package
#   <root>/warp-tools/3rd-party-tools/mmidas/  ← this script
MMIDAS_SRC_DEFAULT="$(cd "$DIR/../../../MMIDAS" 2>/dev/null && pwd)" || true

# Necessary tools and help text
TOOLS=(docker gcloud)
HELP="$(basename "$0") [-h|--help] [-d|--docker-version] [-m|--mmidas-version] [-s|--mmidas-src <path>] [--cache|--no-cache] [-t|--tools] -- build the mmidas image and push to GCR

where:
    -h|--help                  Show help text
    -d|--docker-version <ver>  Image version tag (default: $DOCKER_IMAGE_VERSION)
    -m|--mmidas-version <ver>  MMIDAS package version (default: $MMIDAS_VERSION)
    -s|--mmidas-src <path>     Path to local MMIDAS source tree (default: $MMIDAS_SRC_DEFAULT)
    --cache                    Enable Docker layer caching (default: on)
    --no-cache                 Disable Docker layer caching (force clean rebuild)
    -t|--tools                 Show tools needed to run script

Prerequisites (first-time setup):
    gcloud auth login
    gcloud auth configure-docker us.gcr.io
    "

function main(){
    for t in "${TOOLS[@]}"; do which "$t" >/dev/null || ok=no; done
    if [[ $ok == no ]]; then
        echo "Missing one of the following tools: "
        for t in "${TOOLS[@]}"; do echo "$t"; done
        exit 1
    fi

    local MMIDAS_SRC="$MMIDAS_SRC_DEFAULT"

    while [[ $# -gt 0 ]]; do
        key="$1"
        case $key in
            -d|--docker-version)
            DOCKER_IMAGE_VERSION="$2"
            shift; shift
            ;;
            -m|--mmidas-version)
            MMIDAS_VERSION="$2"
            shift; shift
            ;;
            -s|--mmidas-src)
            MMIDAS_SRC="$2"
            shift; shift
            ;;
            --cache)
            CACHING="ON"
            shift
            ;;
            --no-cache)
            CACHING="OFF"
            shift
            ;;
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

    if [[ -z "$MMIDAS_SRC" || ! -d "$MMIDAS_SRC" ]]; then
        echo "ERROR: MMIDAS source directory not found: '${MMIDAS_SRC}'"
        echo "       Pass -s/--mmidas-src <path> to specify its location."
        exit 1
    fi

    # Check gcloud Docker credential helper is configured for GCR
    if ! gcloud auth print-access-token &>/dev/null; then
        echo "ERROR: gcloud is not authenticated. Run: gcloud auth login"
        exit 1
    fi
    if ! docker-credential-gcloud version &>/dev/null 2>&1; then
        echo "WARNING: GCR credential helper not detected."
        echo "         If the push fails, run: gcloud auth configure-docker us.gcr.io"
    fi

    # ------------------------------------------------------------------
    # Stage the local MMIDAS source into the build context so that
    # Dockerfile can COPY it.  Cleaned up via trap regardless of outcome.
    # ------------------------------------------------------------------
    STAGED="$DIR/mmidas_src"
    trap 'rm -rf "$STAGED"' EXIT

    echo "Staging MMIDAS source from: $MMIDAS_SRC"
    cp -r "$MMIDAS_SRC" "$STAGED"

    IMAGE_TAG="$DOCKER_IMAGE_VERSION-$MMIDAS_VERSION-$TIMESTAMP"
    FULL_IMAGE="$GCR_URL:$IMAGE_TAG"

    echo "Building GCR image: $FULL_IMAGE"
    if [[ "$CACHING" == "ON" ]]; then
        echo "(layer caching enabled)"
        docker build \
            --platform linux/amd64 \
            --build-arg MMIDAS_VERSION="$MMIDAS_VERSION" \
            -t "$FULL_IMAGE" \
            "$DIR"
    else
        docker build --no-cache \
            --platform linux/amd64 \
            --build-arg MMIDAS_VERSION="$MMIDAS_VERSION" \
            -t "$FULL_IMAGE" \
            "$DIR"
    fi

    echo "Pushing to GCR..."
    docker push "$FULL_IMAGE"

    echo -e "$FULL_IMAGE" >> "$DIR/docker_versions.tsv"
    echo "Done: $FULL_IMAGE"
}

main "$@"
