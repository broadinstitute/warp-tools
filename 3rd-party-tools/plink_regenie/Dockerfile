FROM ubuntu:latest 

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" 

WORKDIR /usr/gitc

# need to apt-get everything for python since we can only copy pip installed packages
RUN apt-get update && apt-get install -yq --no-install-recommends --reinstall ca-certificates \
  python-tk \
  tk-dev \
  libssl-dev \
  xz-utils \
  libhdf5-dev \
  openssl \
  make \
  curl \
  unzip \
  libgomp1 \
  liblzo2-dev \
  zlib1g-dev \
  libz-dev \
  git-lfs \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

# Plink install 
RUN mkdir -p /tmp/plink2 && \
  cd /tmp/plink2 && \
  curl -fsSL -o plink2.zip "https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip" && \
  unzip plink2.zip && \
  mv plink2 /bin/plink2 && \
  cd $HOME && \
  rm -rf /tmp/plink2

# REGENIE install
RUN mkdir -p /tmp/regenie && \
  cd /tmp/regenie && \
  curl -fsSL -o regenie.zip "https://github.com/rgcgithub/regenie/releases/download/v2.0.2/regenie_v2.0.2.gz_x86_64_Linux.zip" && \
  unzip regenie.zip && \
  mv regenie_v2.0.2.gz_x86_64_Linux /bin/regenie && \
  cd $HOME && \
  rm -rf /tmp/regenie

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
