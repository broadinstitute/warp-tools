# Image for all Hail-based mtSwirl steps: sharding, merging, combining VCFs,
# finalizing MT, and adding annotations.
#
# Replaces:
#   - warp-tools/3rd-party-tools/aou_mitochondrial_annotate_coverage/Dockerfile (Image 1)
#   - warp/all_of_us/mitochondria/mtSwirl_refactor/Terra/Dockerfile (Image 3)
#
# Build from this directory:
#   docker build -t us.gcr.io/broad-gotc-prod/aou-mito-hail-processing:1.0.0 .

FROM hailgenetics/hail:0.2.127-py3.11

WORKDIR /opt

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libpq-dev \
    git \
    ca-certificates \
    curl \
    gnupg \
    apt-transport-https \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip and handle psycopg2 dependencies
RUN pip install --no-cache-dir pip==25.2 && \
    pip uninstall -y psycopg2 psycopg2-binary || true && \
    pip install psycopg2-binary

# Install Python dependencies
RUN pip install \
    pipdeptree \
    slack-sdk \
    slackclient \
    requests \
    google-cloud-storage \
    --no-cache-dir aiohttp==3.12.15 \
    --no-cache-dir gnomad==0.8.2

# Verify slack installation
RUN python3 -c "from slack import WebClient; print('Slack SDK installed successfully')"

# Verify gnomad_mitochondria installation
RUN python3 -c "import sys; sys.path.append('/opt/gnomad-mitochondria'); import gnomad_mitochondria; print('gnomad_mitochondria available')" || \
    echo "gnomad_mitochondria not installable as package, but should be accessible via PYTHONPATH"

# Clone and install gnomad_qc
RUN git clone https://github.com/broadinstitute/gnomad_qc.git && \
    cd gnomad_qc && \
    pip install . && \
    cd /opt

# Clone and install gnomad-mitochondria
RUN git clone https://github.com/rahulg603/gnomad-mitochondria.git && \
    cd gnomad-mitochondria && \
    ls -la && \
    if [ -f "setup.py" ]; then \
        echo "Found setup.py, installing..."; \
        pip install .; \
    elif [ -f "pyproject.toml" ]; then \
        echo "Found pyproject.toml, installing..."; \
        pip install .; \
    else \
        echo "No setup.py or pyproject.toml found, adding to PYTHONPATH"; \
        echo "export PYTHONPATH=/opt/gnomad-mitochondria:\$PYTHONPATH" >> ~/.bashrc; \
    fi && \
    cd /opt

# Add Google Cloud SDK
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
        > /etc/apt/sources.list.d/google-cloud-sdk.list \
    && curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg \
        | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg \
    && apt-get update \
    && apt-get install -y --no-install-recommends google-cloud-cli \
    && rm -rf /var/lib/apt/lists/*

# Install h5py and gcsfs for coverage DB access
RUN python3 -m pip install --no-cache-dir \
    h5py \
    gcsfs

# Copy mtSwirl source (replaces the git clone of mtSwirl from Image 1)
COPY scripts/ /opt/mtSwirl/generate_mtdna_call_mt/

ENV PYTHONPATH=/opt/mtSwirl

ENTRYPOINT ["python3"]
