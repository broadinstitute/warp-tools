FROM hailgenetics/hail:0.2.133

# Switch to root user to ensure we have permissions
USER root

# Set working directory
WORKDIR /opt

# Upgrade pip and install required Python packages
# Install psycopg2-binary first to avoid compilation issues
RUN pip install --upgrade pip && \
    pip install --upgrade aiohttp && \
    pip install psycopg2-binary && \
    pip install --no-deps gnomad || \
    (echo "Trying alternative installation..." && \
     pip install --upgrade gnomad --force-reinstall --no-cache-dir)

# Clone and setup mtSwirl
RUN git clone https://github.com/jamesemery/mtSwirl.git && \
    chmod -R 755 mtSwirl

# Clone and setup gnomad_qc
RUN git clone https://github.com/broadinstitute/gnomad_qc.git && \
    cd gnomad_qc && \
    pip install . && \
    cd ..

# Clone gnomad-mitochondria (keeping the structure from your original command)
RUN git clone https://github.com/rahulg603/gnomad-mitochondria.git && \
    mv gnomad-mitochondria gnomad_mitochondria_hold && \
    mv gnomad_mitochondria_hold/gnomad_mitochondria ./ && \
    rm -rf gnomad_mitochondria_hold

# Set the working directory back to root for runtime
WORKDIR /

# Add mtSwirl to PATH for easier access
ENV PATH="/opt/mtSwirl:${PATH}"