#!/bin/bash

wget https://github.com/khajoue2/libStatGen/archive/refs/tags/v1.0.15.broad.tar.gz && \
tar -zxvf v1.0.15.broad.tar.gz && \
mv libStatGen-1.0.15.broad libStatGen && \
make -C libStatGen && \

wget http://www.cs.unc.edu/Research/compgeom/gzstream/gzstream.tgz && \
tar -xvf gzstream.tgz && \
make -C gzstream && \
echo "" && \
echo "" && \
echo "" && \
echo "If you are reading this, then the preparations succeeded." && \
echo "You should now be able to run make in this directory."

mkdir -p bin obj
