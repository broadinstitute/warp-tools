# Adding a platform tag to ensure that images built on ARM-based machines doesn't break pipelines
FROM --platform="linux/amd64" python:3.10.12-alpine3.18

ARG STAR_VERSION=2.7.11a

ENV TERM=xterm-256color \
        STAR_URL=https://github.com/alexdobin/STAR/archive/refs/tags/${STAR_VERSION}.tar.gz

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org" \
        STAR_VERSION=${STAR_VERSION}

RUN set -eux; \
        apk add --no-cache \
            bash \
            tini \
            wget \
    ; \
# Install STAR
    mkdir temp; \
    cd temp; \
    \
    wget ${STAR_URL}; \
    tar -xzf ${STAR_VERSION}.tar.gz; \
    cd STAR-${STAR_VERSION}/bin/Linux_x86_64_static; \
    \
    chmod +x ./STAR; \
    cp ./STAR /usr/local/bin; \
    mkdir -p /scripts && \
    apt-get update && apt-get upgrade -y && apt-get install -y libhdf5-dev vim apt-utils liblzma-dev libbz2-dev tini && \
    pip install --upgrade pip && \
    pip install loompy==3.0.6 anndata==0.7.8 numpy==1.23.0 pandas==1.3.5 scipy pysam==0.21

COPY . /scripts
WORKDIR /warptools
# Set tini as default entrypoint
ENTRYPOINT ["/sbin/tini", "--"]
