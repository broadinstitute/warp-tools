FROM python:3.9-slim

# Install system dependencies and tini
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl git libpq-dev build-essential \
    && rm -rf /var/lib/apt/lists/*

ENV TINI_VERSION v0.19.0

RUN wget -O /sbin/tini https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini \
    && chmod +x /sbin/tini

# Pin Python packages: pip, psycopg2-binary, pipdeptree, aiohttp, gnomad
RUN pip install --no-cache-dir pip==25.2

# Clone mtSwirl repo
WORKDIR /usr/gitc
RUN rm -rf mtSwirl \
    && git clone https://github.com/jamesemery/mtSwirl.git

# Clone gnomad_qc, move, and install
RUN rm -rf gnomad_qc gnomad_qc_hold \
    && git clone https://github.com/broadinstitute/gnomad_qc.git \
    && mv gnomad_qc gnomad_qc_hold \
    && mv gnomad_qc_hold/gnomad_qc ./ \
    && rm -rf gnomad_qc_hold \
    && pip install --no-cache-dir ./gnomad_qc

# Clone gnomad-mitochondria, move
RUN rm -rf gnomad_mitochondria gnomad_mitochondria_hold \
    && git clone https://github.com/rahulg603/gnomad-mitochondria.git \
    && mv gnomad-mitochondria gnomad_mitochondria_hold \
    && mv gnomad_mitochondria_hold/gnomad_mitochondria ./ \
    && rm -rf gnomad_mitochondria_hold

# Pip installs & uninstall/reinstall psycopg2 packages to avoid build issues
RUN pip install --no-cache-dir psycopg2-binary==2.9.10 pipdeptree==2.28.0 aiohttp==3.12.15 gnomad==0.8.2 \
    && pip uninstall -y psycopg2 psycopg2-binary \
    && pip install --no-cache-dir psycopg2-binary==2.9.10

# Verify psycopg2 deps via pipdeptree (optional, harmless if it returns nothing)
RUN pip install --no-cache-dir pipdeptree==2.28.0 \
    && pipdeptree | grep psycopg2 || true

ENTRYPOINT ["/sbin/tini", "--"]
WORKDIR /usr/gitc

