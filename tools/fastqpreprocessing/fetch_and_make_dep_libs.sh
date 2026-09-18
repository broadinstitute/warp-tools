#!/bin/bash

wget https://github.com/khajoue2/libStatGen/archive/refs/tags/v1.0.15.broad.tar.gz && \
tar -zxvf v1.0.15.broad.tar.gz && \
mv libStatGen-1.0.15.broad libStatGen && \
make -C libStatGen && \

wget https://github.com/kanedo/gzstream/archive/9a20658673492c41fd9726d439cdea7f96235e84.tar.gz -O gzstream.tar.gz && \
tar -xvf gzstream.tar.gz && \
mv gzstream-9a20658673492c41fd9726d439cdea7f96235e84 gzstream && \
make -C gzstream && \
echo "" && \
echo "" && \
echo "" && \
echo "If you are reading this, then the preparations succeeded." && \
echo "You should now be able to run make in this directory."

mkdir -p bin obj
