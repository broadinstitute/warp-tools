# Adding a platform tag to ensure that images built on ARM-based machines doesn't break pipelines
FROM --platform="linux/amd64" ubuntu:20.04

ARG GLIMPSE_COMMIT_HASH=2a1a8953ac04c512b79cf5ec523045922894ab73

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" \
        GLIMPSE_COMMIT_HASH=${GLIMPSE_COMMIT_HASH}

WORKDIR /docker_build/

# Install required packages
RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -y \
        build-essential \
        libbz2-dev \
        libcurl4-openssl-dev \
        autoconf \
        libssl-dev \
        wget \
        zlib1g-dev \
        liblzma-dev \
        libdeflate-dev \
        git \
        openjdk-8-jre-headless \
        curl \
        python3

# Download and build boost program_options and iostreams
RUN wget https://archives.boost.io/release/1.78.0/source/boost_1_78_0.tar.gz && \
  tar -xf boost_1_78_0.tar.gz && \
  rm boost_1_78_0.tar.gz && \
  cd boost_1_78_0/ && \
  ./bootstrap.sh --with-libraries=iostreams,program_options,serialization --prefix=../boost && \
  ./b2 install && \
  cd .. && \
  cp boost/lib/libboost_iostreams.a boost/lib/libboost_program_options.a boost/lib/libboost_serialization.a /usr/local/lib/ && \
  cp -r boost/include/boost/ /usr/include/ && \
  rm -r boost_1_78_0 boost

# Download and build htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 && \
  tar -xf htslib-1.16.tar.bz2 && \
  rm htslib-1.16.tar.bz2 && \
  cd htslib-1.16 && \
  autoheader && \
  autoconf && \
  ./configure --enable-libcurl && \
  make install && \
  cd .. && \
  rm -r htslib-1.16

# Download and build bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2 && \
  tar -xf bcftools-1.16.tar.bz2 && \
  rm bcftools-1.16.tar.bz2 && \
  cd bcftools-1.16 && \
  autoheader && \
  autoconf && \
  ./configure --enable-libcurl && \
  make install && \
  cd .. && \
  rm -r bcftools-1.16

# Clone GLIMPSE at specific commit hash and build
RUN git clone --recursive https://github.com/odelaneau/GLIMPSE.git && \
    cd GLIMPSE && \
    git checkout ${GLIMPSE_COMMIT_HASH} && \
    make clean && \
    make system -j$(nproc) && \
    mv chunk/bin/GLIMPSE2_chunk \
       split_reference/bin/GLIMPSE2_split_reference \
       phase/bin/GLIMPSE2_phase \
       ligate/bin/GLIMPSE2_ligate \
       concordance/bin/GLIMPSE2_concordance \
       /bin && \
    chmod +x /bin/GLIMPSE2* && \
    cd .. && \
    rm -rf GLIMPSE

# Download picard (for updating sequence dictionary)
RUN wget https://github.com/broadinstitute/picard/releases/download/2.26.11/picard.jar && \
    mv picard.jar /

# Install gcloud
RUN curl -O https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-linux-x86_64.tar.gz && \
    tar -xzf google-cloud-cli-linux-x86_64.tar.gz -C /root && \
    /root/google-cloud-sdk/install.sh --quiet --path-update false --usage-reporting false && \
    rm google-cloud-cli-linux-x86_64.tar.gz

# Add gcloud to PATH
ENV PATH=$PATH:/root/google-cloud-sdk/bin

# Clean up build directory
RUN rm -r /docker_build/


