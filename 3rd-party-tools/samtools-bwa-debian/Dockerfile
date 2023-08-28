# Adding a platform tag to ensure that images built on ARM-based machines doesn't break pipelines
FROM --platform="linux/amd64" debian-slim

ARG SAMTOOLS_VERSION=1.11

ENV TERM=xterm-256color \
    SAMTOOLS_URL=https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        SAMTOOLS_VERSION=${SAMTOOLS_VERSION}

WORKDIR /usr

# Install dependencies
RUN set -eux; \
    apt update; \
    apt-get install -y autoconf automake bash bzip2 gcc wget make libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev numactl; \
    apt-get install -y build-essential git; \
    
    # Install samtools
    mkdir temp; \
    cd temp; \
    wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2; \
    tar -xf samtools-1.11.tar.bz2; \
    cd samtools-1.11; \
    ./configure; \
    make; \
    make install; \
    
    # Install bwa
    cd /usr/temp; \
    git clone --recursive https://github.com/bwa-mem2/bwa-mem2; \
    cd bwa-mem2; \ 
    make; \
    cp bwa-mem2.sse42 /usr/local/bin;\
    cd /usr/local/bin; \
    mv bwa-mem2.sse42 bwa-mem2; \

    # Install tini
    wget https://github.com/krallin/tini/releases/download/v0.19.0/tini -O/sbin/tini; \
    chmod +x /sbin/tini; \
   
    # Clean up cached files
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
