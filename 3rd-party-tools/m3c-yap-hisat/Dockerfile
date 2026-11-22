# Adding a platform tag to ensure that images built on ARM-based machines doesn't break pipelines.
FROM --platform="linux/amd64" mambaorg/micromamba:0.23.0

WORKDIR /usr/gitc
RUN echo $PATH
ENV PATH="$PATH:/opt/conda/bin:/root/pkg/hisat-3n"
RUN echo $PATH

USER root

# Debian 11 (bullseye) is EOL and its security pool has been pruned, so the
# bullseye-security index advertises .deb versions that return 404. Point apt
# at archive.debian.org, which keeps the frozen pool for EOL releases.
# Without this the image no longer builds at all.
RUN printf 'deb http://archive.debian.org/debian bullseye main\ndeb http://archive.debian.org/debian-security bullseye-security main\n' \
        > /etc/apt/sources.list && \
    echo 'Acquire::Check-Valid-Until "false";' > /etc/apt/apt.conf.d/99no-check-valid

# default argument when not provided in the --build-arg
# to build the image with gcp, use
# docker build --build-arg gcp=true -t mapping-gcp:tag .
ARG gcp
RUN if [ "$gcp" = "true" ] ; then \
        apt-get update && \
        apt-get install -y curl gnupg && \
        echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | \
        tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
        curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | \
        apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
        apt-get update -y && \
        apt-get install google-cloud-sdk -y;  \
      else echo 'no gcp install';  \
    fi

# Was `git-all`, which pulls Emacs and Apache/gitweb; its elpa-magit
# post-install script now fails and breaks the build. Only `git clone` is
# needed below, so plain `git` suffices.
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        make \
        git \
        build-essential
    
COPY --chown=$MAMBA_USER:$MAMBA_USER hisat3n_env.yml /tmp/hisat3n_env.yml
RUN micromamba install -y -f /tmp/hisat3n_env.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# cytozip (.cz): compact columnar replacement for allc -- https://github.com/DingWB/cytozip
# Builds Cython extensions with -fopenmp and links against htslib/libdeflate.
# Do NOT add libcurl4-openssl-dev: libhts-dev depends on libcurl4-gnutls-dev and they conflict.
ARG CYTOZIP_VERSION=0.4.5
RUN apt-get update && apt-get install -y --no-install-recommends \
        libhts-dev libdeflate-dev zlib1g-dev libbz2-dev liblzma-dev && \
    rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir cytozip==${CYTOZIP_VERSION}
RUN czip --version && czip bam_to_cz --help > /dev/null

# These run after the cytozip install on purpose: if pip disturbs the pinned
# hisat3n env, the build fails here instead of shipping a broken allcools.
RUN yap --version
RUN allcools --version    
RUN mkdir -p ~/pkg; \
    cd ~/pkg; \
    git clone https://github.com/DaehwanKimLab/hisat2.git hisat-3n; \
    cd hisat-3n; \
    git checkout -b hisat-3n origin/hisat-3n; \
    make; \
    echo 'export PATH=$HOME/pkg/hisat-3n:$PATH' >> ~/.bashrc; \
    source ~/.bashrc
