# Adding a platform tag to ensure that images built on ARM-based machines doesn't break pipelines
FROM --platform="linux/amd64" us.gcr.io/broad-gotc-prod/samtools:2.0.0

ARG BWA_VERSION=2.0pre2

ENV TERM=xterm-256color \
        BWA_URL=https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWA_VERSION}/bwa-mem2-${BWA_VERSION}_x64-linux.tar.bz2
        
LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org" \
        BWA_VERSION=${BWA_VERSION}

WORKDIR /usr/gitc

RUN set -eux; \
        apk add --no-cache \
            gcompat \
            libc6-compat \
            openjdk8 \
            curl \
            git \
            build-base \
    ; \
    # Install bwa-mem2-lisa \
    mkdir temp; \
    cd temp; \
    git clone --recursive https://github.com/bwa-mem2/bwa-mem2 -b bwa-mem2-lisa bwa-mem2-lisa; \
    cd bwa-mem2-lisa; \
    make; \
    cp bwa-mem2.sse42 /usr/local/bin; \
    cd /usr/local/bin; \
    mv bwa-mem2.sse42 bwa-mem2; \
    # Install rust
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh; \
    source "$HOME/.cargo/env";

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
