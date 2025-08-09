FROM python:3.9-slim

# Install system deps
RUN apt-get update && apt-get install -y --no-install-recommends \
    git wget libpq-dev build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install tini
ENV TINI_VERSION=v0.19.0
RUN wget -O /sbin/tini https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini \
    && chmod +x /sbin/tini

WORKDIR /usr/gitc

# Clone mtSwirl repo
RUN rm -rf mtSwirl \
    && git clone https://github.com/jamesemery/mtSwirl.git

# Clone and install gnomad_qc properly (no nested mv, just clone & install)
RUN rm -rf gnomad_qc \
    && git clone https://github.com/broadinstitute/gnomad_qc.git \
    && pip install --no-cache-dir ./gnomad_qc

# Clone gnomad-mitochondria repo
RUN rm -rf gnomad_mitochondria \
    && git clone https://github.com/rahulg603/gnomad-mitochondria.git \
    && mv gnomad-mitochondria gnomad_mitochondria

# Pin pip and python packages with your versions
RUN pip install --no-cache-dir pip==25.2 \
    && pip uninstall -y psycopg2 psycopg2-binary || true \
    && pip install --no-cache-dir psycopg2-binary==2.9.10 pipdeptree==2.28.0 aiohttp==3.12.15 gnomad==0.8.2

# Optional: verify psycopg2 dependency tree (harmless if no output)
RUN pipdeptree | grep psycopg2 || true

# Set tini as default entrypoint
ENTRYPOINT ["/sbin/tini", "--"]

CMD ["python3"]
