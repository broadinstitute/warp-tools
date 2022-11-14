# Adding a platform tag to ensure that images built on ARM-based machines (ex. M-series macs) won't cause issues with our automated PR test suite.
# However, this is not relevant for automated builds in a CI/CD pipeline that is AMD-based.
FROM --platform="linux/amd64" python:3.7.7

# Add metadata labels
LABEL maintainer="Farzaneh Khajouei <fkhajoue@broadinstitute.org>" \
  software="sctools  v.1.0.0" \
  description="A collection of tools for single cell data. Splitting fastq files based on cellbarcodes and other tools to compute metrics on single cell data using barcodes and UMIs."

# Set parameters
ARG htslib_version="1.13"

RUN apt-get update && apt-get upgrade -y && apt-get install -y patch libhdf5-dev vim apt-utils tini &&\
    mkdir /sctools/

COPY ../../../scripts/fastqpreprocessing /sctools/fastqpreprocessing 

RUN cd /sctools/fastqpreprocessing &&\
    wget https://github.com/khajoue2/libStatGen/archive/refs/tags/v1.0.15.broad.tar.gz &&\
    wget https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2 &&\
    tar -zxvf v1.0.15.broad.tar.gz &&\
    tar -jxvf htslib-${htslib_version}.tar.bz2 &&\
    mv libStatGen-1.0.15.broad libStatGen &&\
    cd /sctools/fastqpreprocessing &&\
    wget http://www.cs.unc.edu/Research/compgeom/gzstream/gzstream.tgz &&\
    tar -xvf gzstream.tgz &&\
    cd /sctools/fastqpreprocessing &&\
    make -C libStatGen &&\
    cd /sctools/fastqpreprocessing && make -C htslib-${htslib_version}/ && make -C gzstream &&\
    cd /sctools/fastqpreprocessing && mkdir bin obj && make install &&\
    cp /sctools/fastqpreprocessing/bin/* /usr/local/bin/

WORKDIR usr/local/bin/sctools

# Set tini as default entrypoint
# ENTRYPOINT ["/usr/bin/tini", "--"]
