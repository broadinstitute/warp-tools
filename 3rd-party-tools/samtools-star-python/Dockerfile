FROM --platform="linux/amd64" python:3.10

ARG STAR_VERSION=2.7.11a
ARG SAMTOOLS_VERSION=1.11

ENV TERM=xterm-256color \
    STAR_URL=https://github.com/alexdobin/STAR/archive/refs/tags/${STAR_VERSION}.tar.gz \
    SAMTOOLS_URL=https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
      STAR_VERSION=${STAR_VERSION}

WORKDIR /usr/gitc

RUN set -eux; \
    # Update and install dependencies
    apt-get update && apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
        libhdf5-dev vim apt-utils liblzma-dev libbz2-dev tini \
        autoconf automake bash bzip2 curl gcc libssl-dev make \
        musl-dev ncurses-dev perl wget xz-utils zlib1g-dev && \
    rm -rf /var/lib/apt/lists/* && \
    \
    # Install STAR
    mkdir -p temp && cd temp && \
    wget -O STAR.tar.gz ${STAR_URL} && \
    tar -xzf STAR.tar.gz && \
    cd STAR-${STAR_VERSION}/bin/Linux_x86_64_static && \
    chmod +x STAR && \
    mv STAR /usr/local/bin && \
    \
    # Install SAMtools
    cd /usr/gitc/temp && \
    wget -O samtools.tar.bz2 ${SAMTOOLS_URL} && \
    tar -xf samtools.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure && make && make install && \
    \
    # Cleanup
    cd /usr/gitc && rm -rf temp && \
    mkdir -p /scripts && \
    pip install --upgrade pip && \
    pip install loompy==3.0.6 anndata==0.7.8 numpy==1.23.0 pandas==1.3.5 scipy pysam==0.21

# Copy scripts
COPY . /scripts

# Set working directory
WORKDIR /warptools
