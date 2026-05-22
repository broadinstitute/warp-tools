# Minimal image to run the v2 coverage DB builder (HDF5 + exact summary metrics)
#
# Build from this directory:
#   docker build -t us.gcr.io/broad-gotc-prod/aou-mito-coverage-db:1.0.1 .

FROM python:3.12-slim

# Keep Python output unbuffered for better Cromwell logs
ENV PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PYTHONPATH=/opt/mtSwirl

WORKDIR /opt/mtSwirl

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
       ca-certificates \
  && rm -rf /var/lib/apt/lists/*

RUN python -m pip install --upgrade pip \
  && python -m pip install \
       numpy \
       pandas \
       h5py \
       gcsfs

COPY scripts/ /opt/mtSwirl/generate_mtdna_call_mt/

ENTRYPOINT ["python"]
