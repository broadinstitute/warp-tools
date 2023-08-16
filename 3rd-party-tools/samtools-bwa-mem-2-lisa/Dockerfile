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
    cp /usr/gitc/temp/bwa-mem2-lisa/bwa-mem2.sse42 /usr/local/bin; \
    cd /usr/local/bin; \
    mv bwa-mem2.sse42 bwa-mem2; \
    # Need to cd into directory with bench-fixed-len-e2e-match scripts, etc
    cd /usr/gitc/temp/bwa-mem2-lisa/ext/TAL; \
    # Had to change this for the make command to work
    sed -i '44s/.*/LIBS=-fopenmp -lm -L. -ltal -Lext\/safestringlib\/ -lsafestring -lz/' Makefile; \
    sed -i '1s/.*/#!bin\/sh/' /usr/gitc/temp/bwa-mem2-lisa/ext/TAL/scripts/build-rmi.linear_spline.linear.sh; \   
    make; \ 
    # Install rust
    cd /usr/gitc/temp; \
    curl --proto '=https' --tlsv1.2 https://sh.rustup.rs > rustup.sh && sh rustup.sh -y; \
    source "$HOME/.cargo/env";

# Add bwa-mem2-lisa installation directory to PATH
ENV PATH /usr/gitc/temp/bwa-mem2-lisa/ext/TAL/scripts/:/usr/gitc/temp/bwa-mem2-lisa/ext/TAL/ext/build-rmi/:$PATH


# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
