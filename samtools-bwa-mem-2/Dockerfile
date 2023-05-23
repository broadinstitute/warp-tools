FROM us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616

ARG BWA_VERSION=2.2.1_x64-linux

ENV TERM=xterm-256color \
        BWA_URL=https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-${BWA_VERSION}.tar.bz2
        
LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org" \
        BWA_VERSION=${BWA_VERSION}

WORKDIR /usr/gitc

RUN set -eux; \
        apk add --no-cache \
            gcompat \
            libc6-compat \
            openjdk8 \
    ; \
# Install bwa-mem2 \
    mkdir temp; \
    cd temp; \
    wget -O "bwa.tar.bz2" ${BWA_URL}; \
    tar -xvjf bwa.tar.bz2; \
    cd bwa-mem2-${BWA_VERSION}; \
    make; \
    cp bwa /usr/local/bin \
    ; 

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
