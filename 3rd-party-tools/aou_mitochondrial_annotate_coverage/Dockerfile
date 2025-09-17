FROM gcr.io/broad-dsde-methods/aou-auxiliary/hail_dataproc_wdl:0.2.134

# Set working directory
WORKDIR /opt

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libpq-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip and handle psycopg2 dependencies
RUN pip install --no-cache-dir pip==25.2  && \
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

# Clone and install mtSwirl - updated to point to the main mtswirl repo rather than the 'jamesemery' fork
RUN git clone https://github.com/rahulg603/mtSwirl.git && \
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

# Add gnomad-mitochondria to PYTHONPATH as fallback
ENV PYTHONPATH="/opt/gnomad-mitochondria:${PYTHONPATH}"

# Verify installations
RUN pip list | grep psycopg2 && \
    pipdeptree | grep psycopg2 || true

# Set the working directory back to root for WDL execution
WORKDIR /

# Add mtSwirl to PATH for easier access
ENV PATH="/opt/mtSwirl/generate_mtdna_call_mt:${PATH}"