## Use micromamba as the base image
FROM mambaorg/micromamba:1.5.3

## Set up environment
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
## Create a new environment and install packages
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    conda-forge::r-tidyverse \
    conda-forge::datatable \
    conda-forge::r-optparse \
    conda-forge::r-janitor \
    bioconda::bioconductor-pcatools \
    conda-forge::r-rnomni \
    bioconda::bioconductor-edger \
    bioconda::r-biomartr \
    bioconda::bioconductor-rtracklayer \
    conda-forge::r-patchwork
COPY scripts/*  . 
## Default command to run when the container starts
CMD ["bash"]
