#!/bin/bash
set -e

# Update version when changes to Dockerfile are made
DOCKER_IMAGE_VERSION=1.0.0
TIMESTAMP=$(date +"%s")
DIR=$(cd "$(dirname "$0")" && pwd)

# Registry
GCR_URL="us.gcr.io/broad-gotc-prod/mmidas"

# MMIDAS version (matches setup.py / pyproject.toml at MMIDAS_GIT_REF)
MMIDAS_VERSION="0.1.0"

# Caching: ON by default for faster iterative development.
# Pass --no-cache for a fully clean rebuild when needed.
CACHING="ON"

# ── MMIDAS package source ────────────────────────────────────────────────────
# The mmidas package is installed from a pinned revision of the fork below, not
# from a local directory, so that any published image can be rebuilt from source
# by anyone. This matters because parts of the pipeline live only in that package
# (the pruning loop, K_selection, and the .h5ad loader that supplies gene names
# for KEGG mapping).
#
# MMIDAS_GIT_REF is resolved to an immutable commit SHA before the build, and
# both are recorded as image LABELs and in docker_versions.tsv.
#
# To pick up new MMIDAS changes: commit and push them to the fork, move or add a
# tag, then build with --mmidas-ref <tag>. Passing a full commit SHA also works.
#
# See "Docker image provenance" in pipelines/wdl/mmidas/dashboard.md in
# broadinstitute/warp for what the fork carries on top of upstream.
MMIDAS_GIT_URL="https://github.com/jessicaway/MMIDAS"
MMIDAS_GIT_REF="warp-v2"

# Necessary tools and help text
TOOLS=(docker gcloud git)
HELP="$(basename "$0") [-h|--help] [-d|--docker-version] [-m|--mmidas-version] [-r|--mmidas-ref <tag|sha>] [--cache|--no-cache] [-t|--tools] -- build the mmidas image and push to GCR

where:
    -h|--help                  Show help text
    -d|--docker-version <ver>  Image version tag (default: $DOCKER_IMAGE_VERSION)
    -m|--mmidas-version <ver>  MMIDAS package version (default: $MMIDAS_VERSION)
    -r|--mmidas-ref <tag|sha>  MMIDAS revision to install (default: $MMIDAS_GIT_REF)
                               Resolved to a commit SHA before building.
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
            -r|--mmidas-ref)
            MMIDAS_GIT_REF="$2"
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

    # ------------------------------------------------------------------
    # Resolve MMIDAS_GIT_REF to an immutable commit SHA.
    #
    # A tag can be moved, so the tag is only what we *record* as intent; the
    # SHA is what the build actually fetches. "^{}" peels an annotated tag to
    # the commit it points at (a bare refs/tags/<t> read would give the tag
    # object's own SHA, which is not a valid archive path).
    # ------------------------------------------------------------------
    if [[ "$MMIDAS_GIT_REF" =~ ^[0-9a-f]{40}$ ]]; then
        MMIDAS_SHA="$MMIDAS_GIT_REF"
        echo "MMIDAS ref is already a commit SHA: $MMIDAS_SHA"
    else
        echo "Resolving MMIDAS ref '$MMIDAS_GIT_REF' at $MMIDAS_GIT_URL ..."
        # Query the peeled tag on its own and take it if present. Passing several
        # patterns in one call will not work: git ls-remote sorts its output by
        # refname, so "refs/tags/X" always precedes "refs/tags/X^{}" and a
        # first-match read would return the annotated *tag object* rather than
        # the commit it points at.
        MMIDAS_SHA=$(git ls-remote "${MMIDAS_GIT_URL}.git" \
                        "refs/tags/${MMIDAS_GIT_REF}^{}" 2>/dev/null | cut -f1)
        if [[ -z "$MMIDAS_SHA" ]]; then
            # Lightweight tag or branch: the ref points straight at a commit.
            MMIDAS_SHA=$(git ls-remote "${MMIDAS_GIT_URL}.git" \
                            "refs/tags/${MMIDAS_GIT_REF}" \
                            "refs/heads/${MMIDAS_GIT_REF}" 2>/dev/null \
                         | head -1 | cut -f1)
        fi
        if [[ -z "$MMIDAS_SHA" ]]; then
            echo "ERROR: could not resolve '${MMIDAS_GIT_REF}' in ${MMIDAS_GIT_URL}."
            echo "       Push the tag or branch first, or pass a full 40-character commit SHA."
            exit 1
        fi
    fi
    echo "  MMIDAS revision: ${MMIDAS_GIT_REF} -> ${MMIDAS_SHA}"

    # Fail here rather than mid-build if the archive is not fetchable.
    if ! curl -sfIL -o /dev/null "${MMIDAS_GIT_URL}/archive/${MMIDAS_SHA}.tar.gz"; then
        echo "ERROR: ${MMIDAS_GIT_URL}/archive/${MMIDAS_SHA}.tar.gz is not reachable."
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

    IMAGE_TAG="$DOCKER_IMAGE_VERSION-$MMIDAS_VERSION-$TIMESTAMP"
    FULL_IMAGE="$GCR_URL:$IMAGE_TAG"
    MMIDAS_PROVENANCE="${MMIDAS_GIT_REF}@${MMIDAS_SHA}"

    BUILD_ARGS=(
        --build-arg MMIDAS_VERSION="$MMIDAS_VERSION"
        --build-arg MMIDAS_GIT_URL="$MMIDAS_GIT_URL"
        --build-arg MMIDAS_GIT_REF="$MMIDAS_GIT_REF"
        --build-arg MMIDAS_SHA="$MMIDAS_SHA"
    )

    echo "Building GCR image: $FULL_IMAGE"
    echo "  mmidas from: ${MMIDAS_GIT_URL}/archive/${MMIDAS_SHA}.tar.gz"
    if [[ "$CACHING" == "ON" ]]; then
        echo "(layer caching enabled)"
        docker build \
            --platform linux/amd64 \
            "${BUILD_ARGS[@]}" \
            -t "$FULL_IMAGE" \
            "$DIR"
    else
        docker build --no-cache \
            --platform linux/amd64 \
            "${BUILD_ARGS[@]}" \
            -t "$FULL_IMAGE" \
            "$DIR"
    fi

    echo "Pushing to GCR..."
    docker push "$FULL_IMAGE"

    # Record the MMIDAS revision alongside the tag: the tag alone does not say
    # which model code the image contains.
    echo -e "$FULL_IMAGE\t$MMIDAS_PROVENANCE" >> "$DIR/docker_versions.tsv"
    echo "Done: $FULL_IMAGE"
    echo "  mmidas revision: $MMIDAS_PROVENANCE"
}

main "$@"
