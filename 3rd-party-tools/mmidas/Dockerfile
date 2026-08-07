# Adding a platform tag to ensure that images built on ARM-based machines don't break pipelines
# hadolint ignore=FromPlatformFlagConstDisallowed
FROM --platform="linux/amd64" python:3.10-slim-bookworm

ENV TERM=xterm-256color

# ── Layer 1: system packages ────────────────────────────────────────────────
# Rarely changes. Kept separate so it caches independently.
RUN set -eux; \
    apt-get update; \
    apt-get install -y tini gcc; \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# ── Layer 2: pip bootstrap ───────────────────────────────────────────────────
RUN pip install --upgrade pip

# ── Layer 3: PyTorch ────────────────────────────────────────────────────────
# Isolated in its own layer because it is >2 GB and changes infrequently.
# Only re-downloaded when the version pin below is changed.
RUN pip install torch==2.3.1

# ── Layer 4: remaining scientific dependencies ───────────────────────────────
# Changes occasionally (e.g. scanpy/anndata upgrades).
RUN pip install \
        anndata==0.10.9 \
        numpy==1.26.4 \
        pandas==2.2.2 \
        scipy==1.13.1 \
        scikit-learn==1.5.1 \
        scanpy==1.10.2 \
        seaborn==0.13.2 \
        matplotlib==3.9.1 \
        h5py==3.11.0 \
        hdf5plugin==4.4.0 \
        tensorboard==2.17.0 \
        flatten-json==0.1.14 \
        toml==0.10.2

WORKDIR /usr/local

# ── Metadata ─────────────────────────────────────────────────────────────────
# ARG/LABEL are placed after the heavy installs so that bumping MMIDAS_VERSION
# does not invalidate the torch or deps layers above.
ARG MMIDAS_VERSION=0.1.0
LABEL MAINTAINER="Broad Institute DSDE dsde-engineering@broadinstitute.org" \
      MMIDAS_VERSION=${MMIDAS_VERSION}

# ── Layer 5: MMIDAS package, pinned ──────────────────────────────────────────
# Installed from a specific commit of the fork rather than a copy of a local
# directory, so the image is reproducible from source. MMIDAS_GIT_REF is the
# human-readable tag; MMIDAS_SHA is the immutable commit it resolved to, which
# is what actually gets fetched. docker_build.sh resolves one to the other and
# records both here and in docker_versions.tsv.
#
# A release tarball rather than git+https://: this image has no git binary, and
# adding one to the apt layer would invalidate the ~2 GB torch layer on every
# rebuild.
ARG MMIDAS_GIT_URL=https://github.com/jessicaway/MMIDAS
ARG MMIDAS_GIT_REF=warp-v1
ARG MMIDAS_SHA
RUN test -n "${MMIDAS_SHA}" || { echo "ERROR: MMIDAS_SHA build-arg is required" >&2; exit 1; } \
    && pip install --no-deps "${MMIDAS_GIT_URL}/archive/${MMIDAS_SHA}.tar.gz" \
    && python3 -c "import mmidas, mmidas.utils, mmidas.augmentation; print('mmidas installed:', mmidas.__file__)"

LABEL MMIDAS_GIT_URL=${MMIDAS_GIT_URL} \
      MMIDAS_GIT_REF=${MMIDAS_GIT_REF} \
      MMIDAS_SHA=${MMIDAS_SHA}

# ── Layer 6: pipeline scripts ────────────────────────────────────────────────
# Most frequently edited — last so earlier layers stay cached.
COPY 01_data_prep.py .
COPY test_01_data_prep.py .
COPY 02a_train_augmenter.py .
COPY test_02a_train_augmenter.py .
COPY 02b_train_mixvae.py .
COPY test_02b_train_mixvae.py .
COPY 03a_evaluate.py .
COPY test_03a_evaluate.py .
COPY 03b_classify.py .
COPY test_03b_classify.py .
COPY 03c_traversal_prep.py .
COPY test_03c_traversal_prep.py .
COPY 04_clusterability.py .
COPY test_04_clusterability.py .
COPY 05_state_traversal.py .
COPY test_05_state_traversal.py .

# Set tini as the entrypoint to handle signals correctly
ENTRYPOINT ["/usr/bin/tini", "--"]
