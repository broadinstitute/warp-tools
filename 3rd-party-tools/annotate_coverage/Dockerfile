# Start from the official Hail image (pick your version)
FROM hailgenetics/hail:0.2.133

ENV TINI_VERSION=v0.19.0

# Install system deps needed for gnomad & aiohttp compilation
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    build-essential \
    python3-dev \
    openjdk-11-jdk \
    git \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip, setuptools, wheel first
RUN pip install --no-cache-dir --upgrade pip setuptools wheel

# Install Python dependencies
RUN pip install --no-cache-dir aiohttp==3.12.15 \
    && pip install --no-cache-dir gnomad==0.8.2

# Install mtSwirl
RUN git clone --branch main --depth 1 https://github.com/jamesemery/mtSwirl.git /opt/mtSwirl

# Install tini
RUN wget https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini -O /sbin/tini \
    && chmod +x /sbin/tini

# Set Python path so mtSwirl scripts are discoverable
ENV PYTHONPATH="/opt/mtSwirl:${PYTHONPATH}"

WORKDIR /usr/gitc

# Set tini as default entrypoint
ENTRYPOINT ["/sbin/tini", "--"]

CMD ["python3"]
