# Adding a platform tag to ensure that images built on ARM-based machines doesn't break pipelines
FROM --platform="linux/amd64" eclipse-temurin:17-jdk-alpine

ARG BEAGLE_COMMIT_HASH=d820c4e21ad0ca55548c658f59848038971b9143 \
    BREF3_VERSION=17Dec24.224 \
    BCFTOOLS_VERSION=1.10.2

ENV TERM=xterm-256color \
    TINI_VERSION=v0.19.0

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        BEAGLE_COMMIT_HASH=${BEAGLE_COMMIT_HASH} \
        BREF3_VERSION=${BREF3_VERSION} \
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
            git \
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
    # Clone BEAGLE_COMMIT_HASH from `tmp-sharing/imp-server` repo and generate Beagle jar
    # beagle runs phasing and imputation \
    BEAGLE_COMMIT_SHORT="${BEAGLE_COMMIT_HASH:0:7}"; \
    git clone https://github.com/tmp-sharing/imp-server.git; \
    cd imp-server; \
    git checkout ${BEAGLE_COMMIT_HASH}; \
    cd ..; \
    javac -cp imp-server/src/ imp-server/src/main/Main.java; \
    jar cfe beagle.${BEAGLE_COMMIT_SHORT}.jar main/Main -C imp-server/src/ ./; \
    rm -rf imp-server \
    ; \
    # bref3 converts a reference panel from vcf to the bref3 format that Beagle needs
       curl -L https://faculty.washington.edu/browning/beagle/bref3.${BREF3_VERSION}.jar > bref3.${BREF3_VERSION}.jar \
    ; \
    # Install tini
    wget https://github.com/krallin/tini/releases/download/$TINI_VERSION/tini -O /sbin/tini; \
    chmod +x /sbin/tini;

# Set tini as default entrypoint
ENTRYPOINT ["/sbin/tini", "--" ]
