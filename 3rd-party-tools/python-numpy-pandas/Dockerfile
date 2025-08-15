FROM python:3.9-slim

ENV TINI_VERSION v0.19.0

# Install system dependencies and pinned Python packages
RUN apt-get update && apt-get install -y --no-install-recommends wget \
    && pip install --no-cache-dir pandas==2.2.3 numpy==1.25.2 \
    # Install tini
    && mkdir /tmp/tini-install \
    && cd /tmp/tini-install \
    && wget https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini -O /sbin/tini \
    && chmod +x /sbin/tini \
    && cd / && rm -rf /tmp/tini-install \
    && apt-get purge -y --auto-remove wget \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/gitc

# Set tini as default entrypoint
ENTRYPOINT ["/sbin/tini", "--"]

