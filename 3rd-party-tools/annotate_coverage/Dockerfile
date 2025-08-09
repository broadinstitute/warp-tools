# Start from the official Hail image (pick your version)
FROM hailgenetics/hail:0.2.133

ENV TINI_VERSION=v0.19.0

# Install Tini + dependencies + mtSwirl
RUN apt-get update && apt-get install -y --no-install-recommends wget \
    && pip install --no-cache-dir pip==25.2 \
    && pip install --no-cache-dir --upgrade aiohttp==3.12.15 gnomad==0.8.2 \
    && git clone --branch main --depth 1 https://github.com/jamesemery/mtSwirl.git /opt/mtSwirl \
    && wget https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini -O /sbin/tini \
    && chmod +x /sbin/tini \
    && apt-get purge -y --auto-remove wget \
    && rm -rf /var/lib/apt/lists/*

# Set Python path so mtSwirl scripts are discoverable
ENV PYTHONPATH="/opt/mtSwirl:${PYTHONPATH}"

WORKDIR /usr/gitc

# Set tini as default entrypoint
ENTRYPOINT ["/sbin/tini", "--"]

CMD ["python3"]

