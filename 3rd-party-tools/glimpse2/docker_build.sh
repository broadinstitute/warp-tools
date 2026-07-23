#!/bin/bash
set -e

# This script replicates the palantir-workflows build_base_and_extension_docker.sh approach
# It does a two-stage build:
# 1. Clones GLIMPSE repo and builds the base image using their Dockerfile
# 2. Builds extension image on top using the Dockerfile in this repo (adds gcloud, bcftools, Picard)

# How to run this script:
# 1. Make sure you have Docker and git installed, and Docker is authenticated with GCR if you want to push.
# 2. Run the script with arguments to specify GLIMPSE repo (optional), full commit hash (required), and image tag (optional).
#    For example:
#       ./docker_build.sh -r https://github.com/odelaneau/GLIMPSE.git -c abc123f -t my-custom-tag --record-tag
#    This will build the image using the specified GLIMPSE commit and tag it with 'my-custom-tag'.
#    If you omit the tag, it will auto-generate one based on version, commit hash, and timestamp.
#    If you omit the --record-tag flag, it won't record the tag to docker_versions.tsv (useful for testing or local builds).

# Variables and defaults
# Update version when changes to Dockerfile are made
DOCKER_IMAGE_VERSION=1.2.0
TIMESTAMP=$(date +"%s")
DIR=$(cd $(dirname $0) && pwd)

# GCR registry
GCR_URL="us.gcr.io/broad-gotc-prod/imputation-glimpse2"

# GLIMPSE configuration
GLIMPSE_REPO="https://github.com/odelaneau/GLIMPSE.git"
GLIMPSE_COMMIT_HASH=""  # Required - must be specified by user
CUSTOM_TAG=""  # Optional tag to use instead of auto-generated one


# Help text
HELP="$(basename "$0") [-h|--help] [-r|--repo] -c|--commit [-t|--tag] [--no-push] [--record-tag] [-y|--yes] -- script to build the Imputation GLIMPSE2 image and push to GCR

This script builds a two-stage Docker image:
1. GLIMPSE2 base image (from their repository's Dockerfile)
2. Extension image (adds bcftools, gcloud, Picard on top)

Requirements:
    - docker (required)
    - git (required)
    - Docker must be authenticated with GCR to push (use 'docker login' or 'gcloud auth configure-docker')

where:
    -h|--help           Show help text
    -r|--repo           GLIMPSE repository URL (default: ${GLIMPSE_REPO})
    -c|--commit         GLIMPSE commit hash to use (required)
    -t|--tag            Docker image tag to use (optional, will auto-generate if not provided)
    --no-push           Build locally without pushing to GCR
    --record-tag        Record the image tag to docker_versions.tsv
    -y|--yes            Skip confirmation prompts
    "

# Cleanup function
cleanup() {
    if [ -d "$TEMP_DIR" ]; then
        echo "Cleaning up temporary directory: $TEMP_DIR"
        rm -rf "$TEMP_DIR"
    fi
}
trap cleanup EXIT

function main(){
    PUSH_IMAGE=true
    SKIP_CONFIRM=false
    RECORD_TAG=false

    # Check for required tools
    if ! which docker >/dev/null 2>&1; then
        echo "ERROR: docker is required to run this script"
        exit 1
    fi

    if ! which git >/dev/null 2>&1; then
        echo "ERROR: git is required to run this script"
        exit 1
    fi

    while [[ $# -gt 0 ]]
    do
    key="$1"
    case $key in
        -r|--repo)
        GLIMPSE_REPO="$2"
        shift
        shift
        ;;
        -c|--commit)
        GLIMPSE_COMMIT_HASH="$2"
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

    # Validate required arguments
    if [ -z "$GLIMPSE_COMMIT_HASH" ]; then
        echo "ERROR: GLIMPSE commit hash is required. Use -c or --commit to specify."
        echo ""
        echo "$HELP"
        exit 1
    fi

    # Create temporary directory now that we know we're building
    TEMP_DIR=$(mktemp -d)

    echo "=========================================="
    echo "Building GLIMPSE2 Docker Image (2 stages)"
    echo "=========================================="
    echo "GLIMPSE Repo:   $GLIMPSE_REPO"
    echo "GLIMPSE Commit: $GLIMPSE_COMMIT_HASH"

    # STAGE 1: Clone GLIMPSE repo and build base image
    echo ""
    echo "=========================================="
    echo "STAGE 1: Building Glimpse2 base image..."
    echo "=========================================="

    echo "Cloning GLIMPSE repository..."
    git clone --recursive "$GLIMPSE_REPO" "$TEMP_DIR/GLIMPSE"
    cd "$TEMP_DIR/GLIMPSE"

    # Checkout specified commit
    echo "Checking out Git commit: $GLIMPSE_COMMIT_HASH"
    git checkout "$GLIMPSE_COMMIT_HASH"

    # Get short commit hash for tagging
    GLIMPSE_COMMIT_SHORT=$(echo "$GLIMPSE_COMMIT_HASH" | cut -c1-7)

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
        IMAGE_TAG="$DOCKER_IMAGE_VERSION-$GLIMPSE_COMMIT_SHORT-$TIMESTAMP"
        FINAL_IMAGE_NAME="$GCR_URL:$IMAGE_TAG"
        echo "Auto-generated tag: $IMAGE_TAG"
    fi

    BASE_IMAGE_NAME="temp_glimpse2_base"

    echo "GLIMPSE Commit: $GLIMPSE_COMMIT_HASH (short: $GLIMPSE_COMMIT_SHORT)"
    echo "Base Image:     $BASE_IMAGE_NAME"
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
    echo "Building GLIMPSE2 base image from their Dockerfile..."

    # Check if Dockerfile exists in GLIMPSE repo
    if [ ! -f "Dockerfile" ]; then
        echo "ERROR: Dockerfile not found in GLIMPSE repository!"
        exit 1
    fi

    docker build -t "$BASE_IMAGE_NAME" \
        --platform=linux/amd64 \
        --no-cache \
        .

    echo "GLIMPSE2 base image built successfully: $BASE_IMAGE_NAME"

    # STAGE 2: Build extension image (adds tools on top)
    echo ""
    echo "=========================================="
    echo "STAGE 2: Building extension image (adding bcftools, gcloud, Picard)..."
    echo "=========================================="

    cd "$DIR"
    docker build -t "$FINAL_IMAGE_NAME" \
        --platform=linux/amd64 \
        --no-cache \
        .

    echo "Extension image built successfully: $FINAL_IMAGE_NAME"

    # Clean up temporary base image
    echo ""
    echo "Cleaning up temporary base image..."
    docker rmi "$BASE_IMAGE_NAME" || true

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
    echo "Final image: $FINAL_IMAGE_NAME"
    echo "GLIMPSE Repo: $GLIMPSE_REPO"
    echo "GLIMPSE Commit: $GLIMPSE_COMMIT_HASH"
    echo "=========================================="
}

main "$@"
