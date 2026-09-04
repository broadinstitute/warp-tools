# Adding a platform tag to ensure that images built on ARM-based machines doesn't break pipelines
FROM --platform="linux/amd64" python:3.12-slim-trixie

ARG SCVI_TOOLS_VERSION=1.2
ENV TERM=xterm-256color
LABEL MAINTAINER="Broad Institute DSDE dsde-engineering@broadinstitute.org" \
      SCVI_TOOLS_VERSION=${SCVI_TOOLS_VERSION}

WORKDIR /usr/local

RUN set -eux; \
    apt-get update; \
    apt-get install -y tini wget gcc zlib1g-dev; \
    pip install --upgrade pip; \
    pip install snapatac2==2.7 scvi-tools==${SCVI_TOOLS_VERSION} scanpy numpy google-cloud-storage "scipy<1.14"; \
    pip install --user scikit-misc; \
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.basic.annotation.gff3.gz; \
    apt-get clean && rm -rf /var/lib/apt/lists/*

COPY multiome_label_transfer.py .
COPY gcs_utils.py .
COPY label_transfer_from_preprocessed.py .

# Set tini as the entrypoint
ENTRYPOINT ["/usr/bin/tini", "--"]