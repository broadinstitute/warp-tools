#!/bin/bash

HTSLIBVERSION="1.16"

wget "https://github.com/samtools/htslib/releases/download/$HTSLIBVERSION/htslib-$HTSLIBVERSION.tar.bz2" && \
tar -jxvf "htslib-$HTSLIBVERSION.tar.bz2" && \
mv "htslib-$HTSLIBVERSION" htslib && \
cd htslib && \
./configure --disable-libcurl && \
make && \
cd .. && \

echo "If you are reading this, then the preparations succeeded." && \
echo "You should now be able to run make in this directory."

mkdir -p bin obj
