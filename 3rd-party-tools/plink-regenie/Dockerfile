FROM ubuntu:22.04 

LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org>" 

WORKDIR /usr

# need to apt-get everything for python since we can only copy pip installed packages
RUN set -eux; \
    apt update; \
    apt-get install -y make curl wget unzip libgomp1; \

    # Plink install 
    mkdir -p /tmp/plink2; \
    cd /tmp/plink2; \
    curl -fsSL -o plink2.zip "https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip"; \
    unzip plink2.zip; \
    mv plink2 /bin/plink2; \
    cd $HOME; \
    rm -rf /tmp/plink2; \

   # Regenie 
    mkdir -p /tmp/regenie; \
    cd /tmp/regenie; \
    curl -fsSL -o regenie.zip "https://github.com/rgcgithub/regenie/releases/download/v2.0.2/regenie_v2.0.2.gz_x86_64_Linux.zip"; \
    unzip regenie.zip; \
    mv regenie_v2.0.2.gz_x86_64_Linux /bin/regenie; \
    cd $HOME; \
    rm -rf /tmp/regenie; \

    # Tini
    cd /tmp; \
    wget https://github.com/krallin/tini/releases/download/v0.19.0/tini -O /usr/bin/tini; \
    chmod +x /usr/bin/tini; \

    # Clean up cached files
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Set tini as default entrypoint
ENTRYPOINT [ "/usr/bin/tini", "--" ]
