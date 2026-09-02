FROM python:3.11-slim

ARG MMC_VERSION=1.7.4

RUN apt-get update && \
    apt-get install -y git build-essential wget curl && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir git+https://github.com/AllenInstitute/cell_type_mapper.git@v${MMC_VERSION}

CMD ["/bin/bash"]
