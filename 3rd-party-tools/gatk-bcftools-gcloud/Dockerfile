# Adding a platform tag to ensure that images built on ARM-based machines doesn't break pipelines
FROM --platform="linux/amd64" eclipse-temurin:8-jdk-alpine

ARG BCFTOOLS_VERSION=1.24 \
        GATK4_VERSION=4.2.6.1

ENV TERM=xterm-256color \
    GATK4_URL=https://github.com/broadinstitute/gatk/releases/download/${GATK4_VERSION}/gatk-${GATK4_VERSION}.zip \
    PATH=$PATH:/usr/gitc/gatk4:/root/google-cloud-sdk/bin

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        BCFTOOLS_VERSION=${BCFTOOLS_VERSION} \
        GATK4_VERSION=${GATK4_VERSION}

WORKDIR /usr/gitc

# Install dependencies
RUN set -eux; \
        apk add --no-cache \
            autoconf \
            automake \
            bash \
            bzip2-dev \
            curl \
            curl-dev \
            g++ \
            gcc \
            make \
            musl-dev \
            python3 \
            tini \
            unzip \
            wget \
            xz-dev \
            zlib-dev \
    ; \
# Install GATK
    mkdir gatk4; \
    mkdir temp; \
    cd temp; \
    \
    wget ${GATK4_URL}; \
    unzip gatk-${GATK4_VERSION}.zip; \
    cp -r gatk-${GATK4_VERSION}/* /usr/gitc/gatk4/; \
    chmod +x /usr/gitc/gatk4/gatk; \
    \
    cd ..; \
    rm -r temp \
    ; \
# Install BCFtools
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2; \
    tar xf bcftools-${BCFTOOLS_VERSION}.tar.bz2; \
    cd bcftools-${BCFTOOLS_VERSION}; \
    \
    ./configure; \
    make; \
    make install; \
    \
    cd ../; \
    rm -r bcftools-${BCFTOOLS_VERSION}; \
    rm bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    ; \
# Install HTSlib
    wget https://github.com/samtools/samtools/releases/download/${BCFTOOLS_VERSION}/samtools-${BCFTOOLS_VERSION}.tar.bz2; \
    tar xf samtools-${BCFTOOLS_VERSION}.tar.bz2; \
    cd samtools-${BCFTOOLS_VERSION}/htslib-${BCFTOOLS_VERSION}; \
    \
    ./configure; \
    make; \
    make install; \
    \
    cd ../..; \
    rm -r samtools-${BCFTOOLS_VERSION}; \
    rm samtools-${BCFTOOLS_VERSION}.tar.bz2 \
    ; \
# Install gcloud CLI
    curl -sSL https://sdk.cloud.google.com | bash

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
