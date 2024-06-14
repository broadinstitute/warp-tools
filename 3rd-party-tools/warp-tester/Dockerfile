
# Use the official Ubuntu LTS base image
FROM ubuntu:24.04

LABEL author="Kevin Palis <palis@broadinsitute.org>"

# Install dependencies, tools, and utilities. Then remove the package list to reduce the image size
RUN apt-get update && apt-get install -y \
    curl \
    gnupg \
    lsb-release \
    software-properties-common \
    unzip \
    tree \
    python3 \
    python3-pip \
    kubectl \
    jq \
    netcat \
    && rm -rf /var/lib/apt/lists/*

# Install cromshell
RUN pip3 install cromshell

# Add any additional environment variables here if needed
# ENV VAR_NAME=value

# Set working directory
WORKDIR /workspace


# Set the entrypoint to bash
ENTRYPOINT ["/bin/bash"]
