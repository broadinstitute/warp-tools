# Adding a platform tag to ensure that images built on ARM-based machines doesn't break pipelines
FROM --platform="linux/amd64" debian:buster-slim

ENV TERM=xterm-256color 

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" 

WORKDIR /usr

# Install dependencies
RUN set -eux; \
    apt update; \
    apt-get install -y autoconf automake bash bzip2 gcc wget make libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev numactl; \
    apt-get install -y build-essential git; \

    # Install tini
    wget https://github.com/krallin/tini/releases/download/v0.19.0/tini -O /usr/bin/tini; \
    chmod +x /usr/bin/tini; \
    
    # Install samtools with libdeflate 
    mkdir temp; \
    cd temp; \
    wget https://zenodo.org/record/5731013/files/htstools-1.14_x64-linux.tar.bz2; \
    tar -xf htstools-1.14_x64-linux.tar.bz2; \
    cd htstools-1.14_x64-linux; \
    cp samtools /usr/local/bin;\
    
    # Install bwa
    cd /usr/temp; \
    git clone --recursive https://github.com/bwa-mem2/bwa-mem2; \
    cd bwa-mem2; \ 
    make; \
    cp bwa-mem2.sse42 /usr/local/bin;\
    cd /usr/local/bin; \
    mv bwa-mem2.sse42 bwa-mem2; \
   
    # Clean up cached files
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Set tini as default entrypoint
ENTRYPOINT [ "/usr/bin/tini", "--" ]
