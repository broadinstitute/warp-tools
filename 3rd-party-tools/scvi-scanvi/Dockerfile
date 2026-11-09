# Adding a platform tag to ensure that images built on ARM-based machines doesn't break pipelines
FROM --platform="linux/amd64" python:3.12-slim-trixie

ARG SCVI_TOOLS_VERSION=1.5.1
ENV TERM=xterm-256color
LABEL MAINTAINER="Broad Institute DSDE dsde-engineering@broadinstitute.org" \
      SCVI_TOOLS_VERSION=${SCVI_TOOLS_VERSION}

WORKDIR /usr/local

COPY --from=ghcr.io/astral-sh/uv:0.12.9 /uv /usr/local/bin/uv
# The pinned stack apart from torch; warp-scanvi-training's arm64 image installs the same file.
COPY requirements.txt .

RUN set -eux; \
    grep -qx "scvi-tools==${SCVI_TOOLS_VERSION}" requirements.txt; \
    apt-get update; \
    apt-get install -y tini wget gcc zlib1g-dev; \
    uv pip install --system --no-cache -r requirements.txt torch==2.15.0; \
    uv pip uninstall --system pip; \
    uv pip freeze --system > /usr/local/versions.txt; \
    apt-get clean && rm -rf /var/lib/apt/lists/*

COPY multiome_label_transfer.py label_transfer_from_preprocessed.py preprocess_query.py \
     scp_output.py scanvi_sample.py validate_inputs.py make_reference_only_model.py model_archive.py \
     test_seed_reproducibility.py test_scp_and_query_mapping.py .

# Set tini as the entrypoint
ENTRYPOINT ["/usr/bin/tini", "--"]