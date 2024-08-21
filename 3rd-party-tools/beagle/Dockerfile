# Adding a platform tag to ensure that images built on ARM-based machines doesn't break pipelines
FROM --platform="linux/amd64" adoptopenjdk/openjdk8:alpine-slim

ARG BEAGLE_VERSION=01Mar24.d36 \
    BREF3_VERSION=22Jul22.46e \
    BCFTOOLS_VERSION=1.10.2

ENV TERM=xterm-256color \
    TINI_VERSION=v0.19.0

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        BEAGLE_VERSION=${BEAGLE_VERSION} \
        BCFTOOLS_VERSION=${BCFTOOLS_VERSION}

WORKDIR /usr/gitc

# Install dependencies
RUN set -eux; \
    apk add --no-cache \
            autoconf \
            automake \
            bash \
            bzip2-dev \
            curl \
            g++ \
            gcc \
            gsl-dev \
            make \
            musl-dev \
            perl \
            perl-dev \
            tini \
            wget \
            xz-dev \
            zlib-dev \
    ; \
# Install BCFTools
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
# Download Beagle jars
    # beagle runs phasing and imputation
    curl -L https://faculty.washington.edu/browning/beagle/beagle.${BEAGLE_VERSION}.jar > beagle.${BEAGLE_VERSION}.jar \
    ; \
    # bref3 converts a reference panel from vcf to the bref3 format that Beagle needs
       curl -L https://faculty.washington.edu/browning/beagle/bref3.${BREF3_VERSION}.jar > bref3.${BREF3_VERSION}.jar \
    ; \
# Install tini
    wget https://github.com/krallin/tini/releases/download/$TINI_VERSION/tini -O /sbin/tini; \
    chmod +x /sbin/tini;

# Set tini as default entrypoint
ENTRYPOINT ["/sbin/tini", "--" ]
