# Adding a platform tag to ensure that images built on ARM-based machines
# (ex. M-series macs) won't cause issues with our automated PR test suite.
# However, this is not relevant for automated builds in a CI/CD pipeline that is AMD-based.
FROM --platform="linux/amd64" python:3.7.7

LABEL maintainer="Farzaneh Khajouei <fkhajoue@broadinstitute.org>" \
  software="sctools  v.1.0.0" \
  description="A collection of tools for single cell data. Splitting fastq files based on cellbarcodes and other tools to compute metrics on single cell data using barcodes and UMIs."


RUN apt-get update && apt-get upgrade -y && apt-get install -y libhdf5-dev vim apt-utils liblzma-dev libbz2-dev
RUN mkdir /sctools/

COPY . /sctools

RUN cd /sctools/fastqpreprocessing && ./fetch_and_make_dep_libs.sh

RUN cd /sctools/fastqpreprocessing && make

RUN cp /sctools/fastqpreprocessing/bin/* /usr/local/bin/

WORKDIR usr/local/bin/sctools

# Set tini as default entrypoint
# ENTRYPOINT ["/usr/bin/tini", "--"]
