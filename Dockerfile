# Adding a platform tag to ensure that images built on ARM-based machines
# (ex. M-series macs) won't cause issues with our automated PR test suite.
# However, this is not relevant for automated builds in a CI/CD pipeline that is AMD-based.
FROM --platform="linux/amd64" python:3.7.7

LABEL maintainer="Broad Institute Pipeline Development Team <pipeline-development@broadinstitute.org"  \
  software="warp-tools  v.1.0.0" \
  description="A collection of tools for WARP pipelines."


RUN apt-get update && apt-get upgrade -y && apt-get install -y libhdf5-dev vim apt-utils liblzma-dev libbz2-dev
RUN mkdir -p /sctools

COPY . /sctools

RUN cd /sctools/fastqpreprocessing && ./fetch_and_make_dep_libs.sh && make && cp /sctools/fastqpreprocessing/bin/* /usr/local/bin/

RUN cd /sctools/TagSort && ./fetch_and_make_dep_libs.sh && make && cp /sctools/TagSort/bin/* /usr/local/bin/

WORKDIR usr/local/bin/sctools

# Set tini as default entrypoint
# ENTRYPOINT ["/usr/bin/tini", "--"]
