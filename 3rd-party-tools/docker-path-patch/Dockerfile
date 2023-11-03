# Adding a platform tag to ensure that images built on ARM-based machines
# (ex. M-series macs) won't cause issues with our automated PR test suite.
FROM --platform="linux/amd64" python:3.12

LABEL maintainer="Broad Institute Pipeline Development Team <pipeline-development@broadinstitute.org"  \
  software="Miniwdl Docker Patch by MSFT" \
  description="This program will parse a WDL (using miniwdl) to find instances where the docker runtime parameter is hard-coded as a string URL or specified as a variable with a default in the caller, and replace those instances with variables that are passed up to the root workflow. This allows the user to easily change docker locations without having to modify the WDL."


#create workdir, update package manager, install prerequisite libraries, install python libraries
RUN set -eux && \
    mkdir -p /miniwdl_docker_patch && \
    apt-get update && apt-get upgrade -y && apt-get install -y vim apt-utils tini && \
    pip install --upgrade pip

RUN set -eux && \
  pip install humanize==4.7.0 ipdb==0.13.13 pandas>=1.3.0 miniwdl==1.11.0 rich==13.5.2 

COPY . /miniwdl_docker_patch

WORKDIR /miniwdl_docker_patch

# Set tini as default entry point
ENTRYPOINT ["/usr/bin/tini", "--"]
