# Adding a platform tag to ensure that images built on ARM-based machines doesn't break pipelines.
FROM --platform=linux/amd64 python:3.7.7

LABEL maintainer="DSP Data Centers <monster@broadinstitute.org>"

ARG CUTADAPT_VERSION=4.4

# get cutadapt
RUN pip install --upgrade pip && \
      pip install cutadapt==${CUTADAPT_VERSION}
