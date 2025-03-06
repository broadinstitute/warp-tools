FROM --platform="linux/amd64" us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616

ARG STAR_VERSION=2.7.11a

ENV TERM=xterm-256color \
    STAR_URL=https://github.com/alexdobin/STAR/archive/refs/tags/${STAR_VERSION}.tar.gz

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
      STAR_VERSION=${STAR_VERSION}

WORKDIR /usr/gitc

RUN set -eux; \
    # Update and install necessary dependencies
    apt-get update && apt-get upgrade -y && \
    apt-get install -y libhdf5-dev vim apt-utils liblzma-dev libbz2-dev tini && \
    \
    # Install STAR
    mkdir temp; \
    cd temp; \
    wget ${STAR_URL}; \
    tar -xzf ${STAR_VERSION}.tar.gz; \
    cd STAR-${STAR_VERSION}/bin/Linux_x86_64_static; \
    chmod +x ./STAR; \
    cp ./STAR /usr/local/bin; \
    \
    # Cleanup
    cd /usr/gitc; \
    rm -r temp; \
    \
    # Install Python and required packages
    apt-get install -y python3 python3-pip && \
    pip install --upgrade pip && \
    pip install loompy==3.0.6 anndata==0.7.8 numpy==1.23.0 pandas==1.3.5 scipy pysam==0.21

# Copy scripts
COPY . /scripts

# Set working directory
WORKDIR /warptools

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
