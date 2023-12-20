FROM --platform="linux/amd64" python:3.10.12-buster

LABEL maintainer="Broad Institute Pipeline Development Team <pipeline-development@broadinstitute.org"  \
  software="upstools  2023.03.03" \
  description="Upstools software for Paired-tag."

ENV PATH $PATH:/root/google-cloud-sdk/bin 

RUN set -eux && \
    mkdir -p /upstools && \
    apt-get update && apt-get upgrade -y && apt-get install -y libhdf5-dev vim apt-utils liblzma-dev libbz2-dev tini

COPY . /upstools
#build C/C++ code and add binaries to usr's bin
RUN cd /upstools/scripts/ && sh make.sh && cp /upstools/scripts/* /usr/local/bin/


WORKDIR /upstools

# Set tini as default entry point
ENTRYPOINT ["/usr/bin/tini", "--"]