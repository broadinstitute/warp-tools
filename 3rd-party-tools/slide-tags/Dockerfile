FROM debian:bookworm
ENV LANG C.UTF-8
ENV LANGUAGE C.UTF-8
ENV LC_ALL C.UTF-8

# Install system packages
RUN apt-get update && apt-get install -y \
    wget curl patch gdebi-core gfortran  \
    libssl-dev libxml2-dev               \
    libhdf5-dev hdf5-tools hdf5-helpers  \
    libpng-dev libtiff5-dev libjpeg-dev  \
    libcairo2-dev libfreetype6-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    gnupg

 # Install Google Cloud SDK
 RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
     curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
     apt-get update -y && \
     apt-get install google-cloud-sdk -y

    # Install Julia
RUN wget -P /opt https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.5-linux-x86_64.tar.gz && \
    tar -xzvf /opt/julia-1.10.5-linux-x86_64.tar.gz -C /opt && \
    rm /opt/julia-1.10.5-linux-x86_64.tar.gz && \
    ln -s /opt/julia-1.10.5/bin/julia /usr/local/bin/julia

# Install Julia packages
RUN julia -e 'using Pkg;                 \
              Pkg.add(["CSV",            \
                       "HDF5",           \
                       "FASTX",          \
                       "ArgParse",       \
                       "CodecZlib",      \
                       "IterTools",      \
                       "StatsBase",      \
                       "DataFrames",     \
                       "StringViews",    \
                       "LinearAlgebra",  \
                       "Combinatorics"])'

# Install R
# https://docs.posit.co/resources/install-r.html
RUN export R_VERSION=4.3.3 &&                                                      \
    curl -O https://cdn.rstudio.com/r/debian-12/pkgs/r-${R_VERSION}_1_amd64.deb && \
    echo y | gdebi r-${R_VERSION}_1_amd64.deb &&                                   \
    rm r-${R_VERSION}_1_amd64.deb &&                                               \
    ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R &&                            \
    ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript

# Install R packages
RUN R -e "install.packages(c('tidyverse', 'dplyr', 'tidyr', 'purrr', 'magrittr',  \
                             'future', 'furrr', 'parallelly',                     \
                             'data.table', 'rlist', 'stringr', 'stringi', 'glue', \
                             'ggrastr', 'cowplot', 'gridExtra',        \
                             'scCustomize', 'viridis', 'viridisLite',             \
                             'Seurat', 'SeuratObject',                            \
                             'rdist', 'dbscan', 'RANN',                           \
                             'jsonlite', 'hdf5r', 'qpdf', 'qs', 'qs2',            \
                             'devtools', 'remotes', 'R.utils', 'optparse'),       \
                             repos='http://cloud.r-project.org')"

# Pin ggplot2 to older version
RUN R -e "remotes::install_version("ggplot2", version = "3.5.2", repos = "http://cloud.r-project.org")"

# Install Bioconductor packages
RUN R -e "if (!require('BiocManager', quietly=TRUE)) {install.packages('BiocManager', repos='http://cloud.r-project.org')}; \
          BiocManager::install(c('rhdf5', 'IRanges', 'SummarizedExperiment'))"
# Install other R packages
RUN R -e "remotes::install_github('broadinstitute/DropSift')"
RUN R -e "devtools::install_github('immunogenomics/presto')"

# Download scripts
#RUN wget https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/slide-tags/spatial-count.jl  \
#         https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/slide-tags/run-positioning.R \
#         https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/slide-tags/positioning.R     \
#         https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/slide-tags/helpers.R         \
#         https://raw.githubusercontent.com/MacoskoLab/Macosko-Pipelines/refs/heads/main/slide-tags/plots.R

ENTRYPOINT ["/bin/bash", "-lc"]
CMD ["/bin/bash", "-i"]