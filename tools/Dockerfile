# Adding a platform tag to ensure that images built on ARM-based machines
# (ex. M-series macs) won't cause issues with our automated PR test suite.
# However, this is not relevant for automated builds in a CI/CD pipeline that is AMD-based.
FROM --platform="linux/amd64" python:3.10.12-buster

LABEL maintainer="Broad Institute Pipeline Development Team <pipeline-development@broadinstitute.org"  \
  software="warp-tools  v.1.0.5" \
  description="A collection of tools for WARP pipelines."

#we use gcloud storage copy to download files from google cloud storage -- will be part of warp-tools-c
#used in fastqpreprocessing, ex:
#gcloud storage cp $read1_fastq_files /cromwell_root/input_fastq
ENV PATH $PATH:/root/google-cloud-sdk/bin 

#separating run commands to utilize some caching - each run takes a significant amount of time when not utilizing caching
#create workdir, update package manager, install prerequisite libraries, install python libraries
RUN set -eux && \
    mkdir -p /warptools && \
    apt-get update && apt-get upgrade -y && apt-get install -y libhdf5-dev vim apt-utils liblzma-dev libbz2-dev tini && \
    pip install --upgrade pip && \
    pip install loompy==3.0.6 anndata==0.7.8 numpy==1.23.0 pandas==1.3.5 scipy pysam==0.21 && \
    curl -sSL https://sdk.cloud.google.com | bash 

COPY . /warptools
#build C/C++ code and add binaries to usr's bin
RUN cd /warptools/fastqpreprocessing && ./fetch_and_make_dep_libs.sh && make && cp /warptools/fastqpreprocessing/bin/* /usr/local/bin/
RUN cd /warptools/TagSort && ./fetch_and_make_dep_libs.sh && make && cp /warptools/TagSort/bin/* /usr/local/bin/

WORKDIR /warptools

# Set tini as default entry point
ENTRYPOINT ["/usr/bin/tini", "--"]
