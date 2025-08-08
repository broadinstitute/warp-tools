# Ensure platform consistency (esp. for ARM-based Macs)
FROM --platform=linux/amd64 debian:bullseye-slim

ENV TERM=xterm-256color
LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>"

WORKDIR /usr

# Install dependencies
RUN set -eux; \
    apt update; \
    apt-get install -y --no-install-recommends \
        autoconf automake bash bzip2 gcc wget make \
        libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev \
        numactl build-essential git unzip \
        python3 python3-pip python3-dev \
        python2.7 python2.7-dev \
        openjdk-11-jre-headless \
        ncbi-blast+; \
    \
    # Install pip for python2.7
    wget https://bootstrap.pypa.io/pip/2.7/get-pip.py; \
    python2.7 get-pip.py; \
    rm get-pip.py; \
    \
    # Install tini
    wget https://github.com/krallin/tini/releases/download/v0.19.0/tini -O /usr/bin/tini; \
    chmod +x /usr/bin/tini; \
    \
    # Clean up
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install libidn11 for MitoFinder blast
RUN apt-get update && apt-get install -y libidn11

# Install required Python packages
RUN pip3 install --no-cache-dir biopython requests

# Install Python 2.7 packages for MitoFinder
RUN pip2.7 install --no-cache-dir pillow==6.2.2

# Install MitoFinder
ENV MITOFINDER_VERSION=v1.4

RUN cd /opt && \
    wget https://github.com/RemiAllio/MitoFinder/archive/refs/tags/${MITOFINDER_VERSION}.zip -O mitofinder.zip && \
    unzip mitofinder.zip && \
    mv MitoFinder-${MITOFINDER_VERSION#v} MitoFinder && \
    cd MitoFinder && \
    chmod +x install.sh && \
    ./install.sh && \
    rm /opt/mitofinder.zip


# Add MitoFinder to PATH
ENV PATH="/opt/MitoFinder:${PATH}"

# Create directory for the script
RUN mkdir -p /opt/cleanome/cleanome

# Copy the modified add_mito.py script (Docker version)
COPY add_mito_internal.py /opt/cleanome/cleanome/add_mito.py

# Make add_mito available in PATH
RUN chmod +x /opt/cleanome/cleanome/add_mito.py && \
    ln -s /opt/cleanome/cleanome/add_mito.py /usr/local/bin/add_mito

WORKDIR /data

# Set tini as the entrypoint and run add_mito.py
ENTRYPOINT ["/usr/bin/tini", "--", "python3", "/opt/cleanome/cleanome/add_mito.py"]