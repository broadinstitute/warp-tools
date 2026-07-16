# ===============================================
# STAGE 1: Builder (bcftools + Rust Compilation)
# ===============================================
FROM --platform=linux/amd64 ubuntu:24.04 AS builder

# 1. Install build dependencies & compilers
RUN DEBIAN_FRONTEND=noninteractive apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    clang \
    cmake \
    curl \
    libbz2-dev \
    libcurl4-openssl-dev \
    libdeflate-dev \
    libgsl-dev \
    liblzma-dev \
    libssl-dev \
    wget \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/src

# 2. Download, compile, and stage BCFtools
ENV BCFTOOLS_VERSION=1.24
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar -vxjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && \
    ./configure --prefix=/usr/local --enable-libcurl --with-libdeflate --enable-gsl && \
    make -j$(nproc) && \
    make install DESTDIR=/tmp/bcftools_stage

# 3. Install Git and the Rust Toolchain
RUN DEBIAN_FRONTEND=noninteractive apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y git && \
    rm -rf /var/lib/apt/lists/*

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal
ENV PATH="/root/.cargo/bin:${PATH}"

# 4. Clone the Rust source code from the Git repository
RUN git clone https://github.com/broadinstitute/lrma-sv-imputation-utils.git /usr/src/lrma-sv-imputation-utils
WORKDIR /usr/src/lrma-sv-imputation-utils/resources

# 5. Compile the 3 Rust binaries
RUN cd pop-glimpse2 && cargo build --release
RUN cd paste-vcfs && cargo build --release
RUN cd extract-bubble-PLs && cargo build --release


# ===============================================
# STAGE 2: Runtime Environment
# ===============================================
FROM ubuntu:24.04

# 1. Install ONLY runtime libraries + GCP CLI
RUN DEBIAN_FRONTEND=noninteractive apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg \
    libbz2-1.0 \
    libcurl4 \
    libdeflate0 \
    libgsl27 \
    liblzma5 \
    perl \
    python3 \
    zlib1g \
    && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" > /etc/apt/sources.list.d/google-cloud-sdk.list \
    && apt-get update && \
    apt-get install -y --no-install-recommends google-cloud-cli \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

# 2. Copy the optimized bcftools from the builder stage
COPY --from=builder /tmp/bcftools_stage/usr/local /usr/local/
ENV BCFTOOLS_PLUGINS="/usr/local/libexec/bcftools"

# 3. Copy the 3 compiled Rust binaries from the builder stage into the system PATH
COPY --from=builder /usr/src/lrma-sv-imputation-utils/resources/pop-glimpse2/target/release/pop-glimpse2 /usr/local/bin/
COPY --from=builder /usr/src/lrma-sv-imputation-utils/resources/paste-vcfs/target/release/paste-vcfs /usr/local/bin/
COPY --from=builder /usr/src/lrma-sv-imputation-utils/resources/extract-bubble-PLs/target/release/extract-bubble-PLs /usr/local/bin/

ENV HOME=/root

# Set the default entrypoint
CMD ["bash"]
