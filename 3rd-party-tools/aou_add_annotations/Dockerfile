FROM hailgenetics/hail:0.2.127-py3.11

# Set working directory
WORKDIR /opt

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libpq-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip and handle psycopg2 dependencies
RUN pip install --upgrade pip && \
    pip uninstall -y psycopg2 psycopg2-binary || true && \
    pip install psycopg2-binary

# Install Python dependencies
RUN pip install \
    pipdeptree \
    --upgrade aiohttp \
    --upgrade gnomad

# Clone and install mtSwirl
RUN git clone https://github.com/jamesemery/mtSwirl.git && \
    chmod +x mtSwirl/generate_mtdna_call_mt/add_annotations.py && \
    ls -la mtSwirl/generate_mtdna_call_mt/add_annotations.py

# Clone and install gnomad_qc
RUN git clone https://github.com/broadinstitute/gnomad_qc.git && \
    cd gnomad_qc && \
    pip install . && \
    cd /opt

# Clone and install gnomad-mitochondria
RUN git clone https://github.com/rahulg603/gnomad-mitochondria.git && \
    cd gnomad-mitochondria && \
    pip install . && \
    cd /opt

# Verify installations
RUN pip list | grep psycopg2 && \
    pipdeptree | grep psycopg2 || true && \
    python3 -c "import gnomad_mitochondria; print('gnomad_mitochondria imported successfully')" && \
    python3 -c "import gnomad_qc; print('gnomad_qc imported successfully')"

# Set the working directory back to root for WDL execution
WORKDIR /

# Add mtSwirl to PATH for easier access
ENV PATH="/opt/mtSwirl/generate_mtdna_call_mt:${PATH}"