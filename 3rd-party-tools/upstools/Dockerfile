FROM --platform="linux/amd64" python:3.10.12-buster

LABEL maintainer="Broad Institute Pipeline Development Team <pipeline-development@broadinstitute.org"  \
  software="upstools  v.1.0.0" \
  description="Upstools software for Paired-tag."

ENV PATH $PATH:/root/google-cloud-sdk/bin 

RUN set -eux && \
    mkdir -p /upstools && \
    apt-get update && apt-get upgrade -y && apt-get install -y libhdf5-dev vim apt-utils liblzma-dev libbz2-dev tini && \
    pip install --upgrade pip && \
    curl -sSL https://sdk.cloud.google.com | bash 

COPY . /scripts
#build C/C++ code and add binaries to usr's bin
RUN cd /scripts/ && ./fetch_and_make_dep_libs.sh && make && cp /upstools/bin/* /usr/local/bin/


WORKDIR /upstools

# Set tini as default entry point
ENTRYPOINT ["/usr/bin/tini", "--"]