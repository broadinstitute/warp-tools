#!/bin/bash
set -e

# This script builds the SV Imputation Rust Tools Docker image and pushes it to GCR.
# The image includes BCFtools and three Rust binaries from the lrma-sv-imputation-utils repository.
# (https://github.com/broadinstitute/lrma-sv-imputation-utils)

# How to run this script:
# 1. Make sure you have Docker installed and Docker is authenticated with GCR if you want to push.
# 2. Run the script with arguments to specify git commit (optional) and image tag (optional).
#    For example:
#       ./docker_build.sh -c abc123def456... -t my-custom-tag --record-tag
#    This will build the image using the specified git commit hash and tag it with 'my-custom-tag'.
#    If you omit the git commit, it will use the latest commit from 'main' branch.
#    If you omit the tag, it will auto-generate one based on version, git commit, and timestamp.
#    If you omit the --record-tag flag, it won't record the tag to docker_versions.tsv (useful for testing or local builds).

# Variables and defaults
# Update version when changes to Dockerfile are made
DOCKER_IMAGE_VERSION=1.0.0
TIMESTAMP=$(date +"%s")
DIR=$(cd $(dirname $0) && pwd)

# GCR registry
#GCR_URL="us.gcr.io/broad-gotc-prod/sv-imputation-rust-tools" # TODO - Saloni uncomment this
GCR_URL="us.gcr.io/broad-dsde-methods/sshah/sv_rust_tools"

# Git commit configuration (must be a full 40-character commit hash, or empty to use latest main)
GIT_COMMIT=""  # Default to empty (will use latest commit from main branch)
CUSTOM_TAG=""  # Optional tag to use instead of auto-generated one

# Help text
HELP="$(basename "$0") [-h|--help] [-c|--commit] [-t|--tag] [--no-push] [--record-tag] [-y|--yes] -- script to build the SV Imputation Rust Tools image and push to GCR

This script builds a Docker image containing:
- BCFtools
- Three Rust binaries from lrma-sv-imputation-utils repository(https://github.com/broadinstitute/lrma-sv-imputation-utils):
  * pop-glimpse2
  * paste-vcfs
  * extract-bubble-PLs

Requirements:
    - docker (required)
    - Docker must be authenticated with GCR to push (use 'docker login' or 'gcloud auth configure-docker')

where:
    -h|--help           Show help text
    -c|--commit         Full 40-character git commit hash for lrma-sv-imputation-utils repo (optional, uses latest main if not provided)
    -t|--tag            Docker image tag to use (optional, will auto-generate if not provided)
    --no-push           Build locally without pushing to GCR
    --record-tag        Record the image tag to docker_versions.tsv
    -y|--yes            Skip confirmation prompts
    "

function main(){
    PUSH_IMAGE=true
    SKIP_CONFIRM=false
    RECORD_TAG=false

    # Check for required tools
    if ! which docker >/dev/null 2>&1; then
        echo "ERROR: docker is required to run this script"
        exit 1
    fi

    while [[ $# -gt 0 ]]
    do
    key="$1"
    case $key in
        -c|--commit)
        GIT_COMMIT="$2"
        # Validate that it's a full 40-character commit hash
        if [ -n "$GIT_COMMIT" ] && ! [[ "$GIT_COMMIT" =~ ^[0-9a-f]{40}$ ]]; then
            echo "ERROR: GIT_COMMIT must be a full 40-character commit hash (lowercase hex)"
            echo "Provided: $GIT_COMMIT"
            exit 1
        fi
        shift
        shift
        ;;
        -t|--tag)
        CUSTOM_TAG="$2"
        shift
        shift
        ;;
        -h|--help)
        echo "$HELP"
        exit 0
        ;;
        --no-push)
        PUSH_IMAGE=false
        shift
        ;;
        --record-tag)
        RECORD_TAG=true
        shift
        ;;
        -y|--yes)
        SKIP_CONFIRM=true
        shift
        ;;
        *)
        echo "ERROR: Unknown argument: $1"
        echo ""
        echo "$HELP"
        exit 1
        ;;
    esac
    done

    echo "=========================================="
    echo "Building SV Imputation Rust Tools Docker Image"
    echo "=========================================="
    if [ -n "$GIT_COMMIT" ]; then
        echo "Git Commit:     $GIT_COMMIT"
    else
        echo "Use latest commit from main branch"
    fi

    # Get short git commit for tagging (use first 7 characters if full hash provided)
    if [ -n "$GIT_COMMIT" ]; then
        GIT_COMMIT_SHORT=$(echo "$GIT_COMMIT" | cut -c1-7)
    else
        # Fetch the latest commit hash from the main branch of the remote repo
        echo "Fetching latest commit hash from main branch..."
        LATEST_COMMIT=$(git ls-remote https://github.com/broadinstitute/lrma-sv-imputation-utils.git refs/heads/main | cut -f1)
        if [ -z "$LATEST_COMMIT" ]; then
            echo "ERROR: Failed to fetch latest commit from main branch"
            exit 1
        fi
        echo "Latest commit on main: $LATEST_COMMIT"
        # Set GIT_COMMIT to the fetched commit so it's passed to Docker
        GIT_COMMIT="$LATEST_COMMIT"
        GIT_COMMIT_SHORT=$(echo "$LATEST_COMMIT" | cut -c1-7)
    fi

    # Use custom tag if provided, otherwise auto-generate
    if [ -n "$CUSTOM_TAG" ]; then
        # Check if custom tag already contains the full image path (registry/path:tag)
        if [[ "$CUSTOM_TAG" == *":"* ]] && [[ "$CUSTOM_TAG" == *"/"* ]]; then
            # Tag already includes registry/path:tag, use it as-is
            FINAL_IMAGE_NAME="$CUSTOM_TAG"
            # Extract just the tag portion for display
            IMAGE_TAG="${CUSTOM_TAG##*:}"
            echo "Using provided full image path: $FINAL_IMAGE_NAME"
        else
            # Tag is just the version/label, prepend GCR URL
            IMAGE_TAG="$CUSTOM_TAG"
            FINAL_IMAGE_NAME="$GCR_URL:$IMAGE_TAG"
            echo "Using provided tag: $IMAGE_TAG"
        fi
    else
        IMAGE_TAG="$DOCKER_IMAGE_VERSION-$GIT_COMMIT_SHORT-$TIMESTAMP"
        FINAL_IMAGE_NAME="$GCR_URL:$IMAGE_TAG"
        echo "Auto-generated tag: $IMAGE_TAG"
    fi

    echo "Registry:       $GCR_URL"
    echo "Final Image:    $FINAL_IMAGE_NAME"
    echo "=========================================="

    # Confirm before building
    if [ "$SKIP_CONFIRM" = false ]; then
        read -p "Continue with build? (y/n) " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "Build cancelled."
            exit 0
        fi
    fi

    echo ""
    echo "Building Docker image..."
    echo ""

    cd "$DIR"
    docker build -t "$FINAL_IMAGE_NAME" \
        --build-arg GIT_COMMIT="$GIT_COMMIT" \
        --platform=linux/amd64 \
        --no-cache \
        .

    echo ""
    echo "Image built successfully: $FINAL_IMAGE_NAME"

    # Push to GCR if requested
    if [ "$PUSH_IMAGE" = true ]; then
        echo ""
        echo "Pushing image to GCR..."
        docker push "$FINAL_IMAGE_NAME"
        echo "Image pushed successfully"
    else
        echo ""
        echo "Skipping push to GCR (--no-push flag set)"
        echo "Image available locally as: $FINAL_IMAGE_NAME"
    fi

    # Record tag if requested
    if [ "$RECORD_TAG" = true ]; then
        echo ""
        echo "$FINAL_IMAGE_NAME" >> "$DIR/docker_versions.tsv"
        echo "Recorded image tag to docker_versions.tsv"
    fi

    echo ""
    echo "=========================================="
    echo "Build complete!"
    echo "Git Commit: $GIT_COMMIT"
    echo "Final image: $FINAL_IMAGE_NAME"
    echo "=========================================="
}

main "$@"
