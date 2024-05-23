# Adding a platform tag to ensure that images built on ARM-based machines doesn't break pipelines
FROM --platform="linux/amd64" python:3.9

ENV TERM=xterm-256color \
        TINI_VERSION=v0.19.0\
        PATH=$PATH:/usr/gitc

ARG SNAPATAC2_VERSION=2.6.3

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org" \
        SNAPATAC2_VERSION=${SNAPATAC2_VERSION}

WORKDIR /usr/gitc

RUN apt-get update; \
    apt-get install -y tabix; \
    pip install cmake; \
    curl https://sh.rustup.rs -sSf | sh -s -- -y; \
    source "$HOME/.cargo/env"; \
    pip install snapatac2==${SNAPATAC2_VERSION}; \
    # Install TINI
    mkdir temp; \
    cd temp; \
    wget https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini -O /sbin/tini; \
    chmod +x /sbin/tini; \
    rm -r /usr/gitc/temp;

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
