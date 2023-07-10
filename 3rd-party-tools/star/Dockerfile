# Adding a platform tag to ensure that images built on ARM-based machines doesn't break pipelines
FROM --platform="linux/amd64" alpine:latest

ARG STAR_VERSION=STAR_2.7.10b_alpha_2023-07-06_Linux_x86_64_static

ENV TERM=xterm-256color \
        STAR_URL=https://github.com/dobinlab/STAR_pre_releases/releases/download/2.7.10b_alpha_2023-07-06/STAR_2.7.10b_alpha_2023-07-06_Linux_x86_64_static

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org" \
        STAR_VERSION=${STAR_VERSION}

WORKDIR /usr/gitc

RUN set -eux; \
        apk add --no-cache \
            bash \
            tini \
            wget \
    ; \
# Install STAR
    mkdir temp; \
    cd temp; \
    \
    wget ${STAR_URL}; \
    tar -xzf ${STAR_VERSION}.tar.gz; \
    cd STAR-${STAR_VERSION}/bin/Linux_x86_64_static; \
    \
    chmod +x ./STAR; \
    cp ./STAR /usr/local/bin; \
    \
    rm -r /usr/gitc/temp;

# Set tini as default entrypoint
ENTRYPOINT ["/sbin/tini", "--"]
