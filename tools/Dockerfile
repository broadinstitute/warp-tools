# Adding a platform tag to ensure that images built on ARM-based machines
# (ex. M-series macs) won't cause issues with our automated PR test suite.
# However, this is not relevant for automated builds in a CI/CD pipeline that is AMD-based.
FROM --platform="linux/amd64" python:3.10.12-buster

LABEL maintainer="Broad Institute Pipeline Development Team <pipeline-development@broadinstitute.org"  \
  software="warp-tools  v.1.0.1" \
  description="A collection of tools for WARP pipelines."

ENV PATH $PATH:/root/google-cloud-sdk/bin 

#separating run commands to utilize some caching - each run takes a significant amount of time when not utilizing caching
#create workdir, update package manager, install prerequisite libraries, install python libraries
RUN set -eux && \
    mkdir -p /warptools && \
    apt-get update && apt-get upgrade -y && apt-get install -y libhdf5-dev vim apt-utils liblzma-dev libbz2-dev && \
    pip install --upgrade pip && \
    pip install loompy==3.0.6 anndata==0.7.8 numpy==1.23.0 pandas==1.3.5 scipy h5py==2.10.0 && \
    curl -sSL https://sdk.cloud.google.com | bash 

COPY . /warptools
#build C/C++ code and add binaries to usr's bin
RUN cd /warptools/fastqpreprocessing && ./fetch_and_make_dep_libs.sh && make && cp /warptools/fastqpreprocessing/bin/* /usr/local/bin/
RUN cd /warptools/TagSort && ./fetch_and_make_dep_libs.sh && make && cp /warptools/TagSort/bin/* /usr/local/bin/

WORKDIR /warptools

# Set tini as default entrypoint
# ENTRYPOINT ["/usr/bin/tini", "--"]
