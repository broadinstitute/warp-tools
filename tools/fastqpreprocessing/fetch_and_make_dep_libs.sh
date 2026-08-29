#!/bin/bash

wget https://github.com/khajoue2/libStatGen/archive/refs/tags/v1.0.15.broad.tar.gz && \
tar -zxvf v1.0.15.broad.tar.gz && \
mv libStatGen-1.0.15.broad libStatGen && \
sed -i 's/-Werror//g' libStatGen/general/Makefile && \
make -C libStatGen && \

wget https://github.com/kanedo/gzstream/archive/refs/heads/master.tar.gz -O gzstream.tar.gz && \
tar -xvf gzstream.tar.gz && \
mv gzstream-master gzstream && \
make -C gzstream && \
echo "" && \
echo "" && \
echo "" && \
echo "If you are reading this, then the preparations succeeded." && \
echo "You should now be able to run make in this directory."

mkdir -p bin obj
