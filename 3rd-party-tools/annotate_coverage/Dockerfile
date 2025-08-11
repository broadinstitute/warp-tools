# Start from the official Hail image
FROM hailgenetics/hail:0.2.133

ENV TINI_VERSION=v0.19.0

# Install system deps + pip
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    build-essential \
    python3-dev \
    python3-pip \
    openjdk-11-jdk \
    git \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip, setuptools, wheel
RUN python3 -m pip install --no-cache-dir --upgrade pip setuptools wheel

# Install Python dependencies
RUN python3 -m pip install --no-cache-dir aiohttp==3.12.15 \
    && python3 -m pip install --no-cache-dir gnomad==0.8.2

# Install mtSwirl
RUN git clone --branch main --depth 1 https://github.com/jamesemery/mtSwirl.git /opt/mtSwirl

# Install tini
RUN wget https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini -O /sbin/tini \
    && chmod +x /sbin/tini

# Set Python path so mtSwirl scripts are discoverable
ENV PYTHONPATH="/opt/mtSwirl:\${PYTHONPATH}"

WORKDIR /usr/gitc

ENTRYPOINT ["/sbin/tini", "--"]
CMD ["python3"]
