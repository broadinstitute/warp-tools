FROM ubuntu:18.04

#### Basic image utilities
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
    python3 \
    wget \
    curl \
    bc \
    unzip \
    bzip2 \
    less \
    bedtools \
    samtools \
    tabix \
    gpg-agent \
    build-essential \
    git \
    openjdk-17-jdk \
    software-properties-common && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove && \
    rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV DOWNLOAD_DIR /downloads
ENV CONDA_URL https://repo.anaconda.com/miniconda/Miniconda3-py310_23.1.0-1-Linux-x86_64.sh 
ENV CONDA_MD5 = "32d73e1bc33fda089d7cd9ef4c1be542616bd8e437d1f77afeeaf7afdb019787"
ENV CONDA_PATH /opt/miniconda
ENV PATH $CONDA_PATH/bin:$PATH
RUN mkdir $DOWNLOAD_DIR && \
    wget -nv -O $DOWNLOAD_DIR/miniconda.sh $CONDA_URL && \
    test "`md5sum $DOWNLOAD_DIR/miniconda.sh | awk -v FS='  ' '{print $1}'` = $CONDA_MD5" && \
    bash $DOWNLOAD_DIR/miniconda.sh -p $CONDA_PATH -b && \
    rm $DOWNLOAD_DIR/miniconda.sh && \
    conda clean -afy && \
    find /opt/miniconda/ -follow -type f -name '*.a' -delete && \
    find /opt/miniconda/ -follow -type f -name '*.pyc' -delete && \
    rm -rf /root/.cache/pip

RUN git clone --depth 1 https://github.com/Illumina/Cyrius.git

RUN pip install -r Cyrius/requirements.txt && conda install -c conda-forge google-cloud-sdk
 

