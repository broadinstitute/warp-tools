# Use the Hail base image (includes hail, numpy, pandas, pyspark, Python 3.11)
# Mirror of hailgenetics/hail to avoid DockerHub pull limits
FROM hailgenetics/hail:0.2.134-py3.11

WORKDIR /opt

# Install additional Python packages needed by mito_plot_filter.py
RUN pip install --no-cache-dir \
    matplotlib \
    seaborn

# Install tabix (required for hl.export_vcf tabix=True)
RUN apt-get update && apt-get install -y --no-install-recommends \
    tabix \
    && rm -rf /var/lib/apt/lists/*

# Copy the post-processing script into the image
COPY mito_plot_filter.py /opt/mito_plot_filter.py

# Smoke-test: confirm all imports resolve
RUN python3 -c "\
import hail as hl; \
import matplotlib; \
import matplotlib.pyplot as plt; \
import numpy as np; \
import pandas as pd; \
import seaborn as sns; \
print('All imports OK')"

CMD ["bash"]
