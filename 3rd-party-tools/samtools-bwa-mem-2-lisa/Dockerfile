# Adding a platform tag to ensure that images built on ARM-based machines doesn't break pipelines
FROM --platform="linux/amd64" us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616

ARG BWA_VERSION=2.0pre2

ENV TERM=xterm-256color \
        BWA_URL=https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWA_VERSION}/bwa-mem2-${BWA_VERSION}_x64-linux.tar.bz2
        
LABEL MAINTAINER="Broad Institute DSDE <dsde-engineering@broadinstitute.org" \
        BWA_VERSION=${BWA_VERSION}

WORKDIR /usr/gitc

RUN set -eux; \
        apk add --no-cache \
            gcompat \
            libc6-compat \
            openjdk8 \
            curl \
            git \
            build-base \
    ; \

    # Install rust
    mkdir temp; \
    cd temp; \
    curl --proto '=https' --tlsv1.2 https://sh.rustup.rs > rustup.sh && sh rustup.sh -y; \
    source "$HOME/.cargo/env";
    
    # Install bwa-mem2-lisa \
    cd /usr/gitc/temp; \
    git clone --recursive https://github.com/bwa-mem2/bwa-mem2 -b bwa-mem2-lisa bwa-mem2-lisa; \
    cd bwa-mem2-lisa; \
    
    # Need to cd into directory with bench-fixed-len-e2e-match scripts, etc
    cd /usr/gitc/temp/bwa-mem2-lisa/ext/TAL; \
    # Had to change this for the make command to work
    sed -i '44s/.*/LIBS=-fopenmp -lm -L. -ltal -Lext\/safestringlib\/ -lsafestring -lz/' Makefile; \
    # Need to update paths in /ext/TAL/scripts/build-rmi.linear_spline.linear.sh and /ext/TAL/src/LISA-FMI/ipbwt_rmi.h
    # replace lines 52, 53, 455, 704
    sed -i '453s#.*#string shell_cmd = "cd " + lisa_home + "&& bash /usr/gitc/temp/bwa-mem2-lisa/ext/TAL/scripts/build-rmi.linear_spline.linear.sh " + ipbwt_f64_filename + " " + rmi_filename + " " + to_string(num_rmi_leaf_nodes) + " " + "F64";#' ipbwt_rmi.h; \
    sed -i '52s#.*#IPBWT_RMI(const string &t, string ref_seq_filename, int K_, int64_t num_rmi_leaf_nodes, index_t *__sa = NULL, string lisa_home = "/usr/gitc/temp/bwa-mem2-lisa");#' ipbwt_rmi.h; \
    sed -i '53s#.*#IPBWT_RMI(const string \&t, index_t t_size, string ref_seq_filename, int K_, int64_t num_rmi_leaf_nodes, index_t *__sa = NULL, string lisa_home = "/usr/gitc/temp/bwa-mem2-lisa");#' ipbwt_rmi.h; \
    
    cd /usr/gitc/temp/bwa-mem2-lisa; 
    make; \
    cp /usr/gitc/temp/bwa-mem2-lisa/bwa-mem2.sse42 /usr/local/bin; \
    cd /usr/local/bin; \
    mv bwa-mem2.sse42 bwa-mem2; \

# Set tini as default entrypoint
ENTRYPOINT [ "/sbin/tini", "--" ]
