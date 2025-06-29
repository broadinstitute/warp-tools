RUN apt-get update; \
    apt-get install -y tabix; \
    pip install cmake; \
    curl https://sh.rustup.rs -sSf | sh -s -- -y; \
    source "$HOME/.cargo/env"; \
    pip install snapatac2==${SNAPATAC2_VERSION}; \
    # Install the CPU-safe version of polars to avoid illegal instruction crash
    pip install "polars-lts-cpu>=0.20.0"; \
    pip install scanpy; \
    # Install TINI
    mkdir temp; \
    cd temp; \
    wget https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini -O /sbin/tini; \
    chmod +x /sbin/tini; \
    rm -r /usr/gitc/temp;
