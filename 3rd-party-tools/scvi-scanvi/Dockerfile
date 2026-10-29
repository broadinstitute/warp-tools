# Adding a platform tag to ensure that images built on ARM-based machines doesn't break pipelines
FROM --platform="linux/amd64" python:3.12-slim-trixie

ARG SCVI_TOOLS_VERSION=1.5.1
ENV TERM=xterm-256color
LABEL MAINTAINER="Broad Institute DSDE dsde-engineering@broadinstitute.org" \
      SCVI_TOOLS_VERSION=${SCVI_TOOLS_VERSION}

WORKDIR /usr/local

RUN set -eux; \
    apt-get update; \
    apt-get install -y tini wget gcc zlib1g-dev; \
    pip install --upgrade pip; \
    pip install snapatac2==2.10.0 scvi-tools==${SCVI_TOOLS_VERSION} scanpy numpy; \
    pip install --user scikit-misc; \
    pip freeze > /usr/local/versions.txt; \
    apt-get clean && rm -rf /var/lib/apt/lists/*

COPY multiome_label_transfer.py label_transfer_from_preprocessed.py preprocess_query.py \
     scp_output.py scanvi_sample.py validate_inputs.py make_reference_only_model.py model_archive.py \
     test_seed_reproducibility.py test_scp_and_query_mapping.py .

# Set tini as the entrypoint
ENTRYPOINT ["/usr/bin/tini", "--"]