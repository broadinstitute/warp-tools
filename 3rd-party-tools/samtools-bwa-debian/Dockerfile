# Adding a platform tag to ensure that images built on ARM-based machines doesn't break pipelines
FROM --platform="linux/amd64" debian

ARG SAMTOOLS_VERSION=1.11

ENV TERM=xterm-256color \
    SAMTOOLS_URL=https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        SAMTOOLS_VERSION=${SAMTOOLS_VERSION}

WORKDIR /usr

# Install dependencies
RUN set -eux; \
    apt update; \
    apt-get install autoconf automake bash bzip2 gcc wget make libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev numactl; \
    apt-get install build-essential git; \

    # Install samtools
    mkdir temp; \
    cd temp; \
    wget ${SAMTOOLS_URL}; \
    tar -xf samtools-${SAMTOOLS_VERSION}.tar.bz2; \
    cd samtools-${SAMTOOLS_VERSION}; \
    \
    ./configure; \
    make; \
    make install; \
    
    # Install bwa
    cd /usr/temp; \
    git clone --recursive https://github.com/bwa-mem2/bwa-mem2; \
    cd bwa-mem2; \ 
    make; \
    sudo cp bwa-mem2.sse42 /usr/local/bin;\
    cd /usr/local/bin; \
    sudo mv bwa-mem2.sse42 bwa-mem2; 


# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
