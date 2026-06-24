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

# ── Layer 5: local MMIDAS package ────────────────────────────────────────────
# Re-runs when the source tree changes (staged by docker_build.sh).
# The upstream setup.py only declares packages=["mmidas"], omitting subpackages.
# Patch it to use find_packages() so mmidas.utils and mmidas.augmentation are
# importable. __init__.py files are present in the fork.
COPY mmidas_src/ /tmp/mmidas_src/
RUN sed -i 's/from setuptools import setup/from setuptools import setup, find_packages/' /tmp/mmidas_src/setup.py \
    && sed -i 's/packages=\["mmidas"\]/packages=find_packages()/' /tmp/mmidas_src/setup.py \
    && pip install --no-deps /tmp/mmidas_src \
    && rm -rf /tmp/mmidas_src

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
