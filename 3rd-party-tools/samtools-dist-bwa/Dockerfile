# Adding a platform tag to ensure that images built on ARM-based machines doesn't break pipelines
FROM --platform="linux/amd64"  ubuntu:22.04

ENV TERM=xterm-256color 

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" 

WORKDIR /usr

# Install dependencies
RUN set -eux; \
    apt update; \
    apt-get install -y autoconf automake bash bzip2 gcc wget make libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev numactl; \
    apt-get install -y build-essential git; \

    # Install open omics acceleration framework
    mkdir temp; \
    ls; \
    pwd; \
    cd temp; \
    git clone --recursive https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git; \
    cd Open-Omics-Acceleration-Framework/pipelines/fq2sortedbam; \
    ./install.sh; \

    # Tini
    cd /temp; \
    wget https://github.com/krallin/tini/releases/download/v0.19.0/tini -O /usr/bin/tini; \
    chmod +x /usr/bin/tini; \
   
    # Clean up cached files
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Set tini as default entrypoint
ENTRYPOINT [ "/usr/bin/tini", "--" ]
