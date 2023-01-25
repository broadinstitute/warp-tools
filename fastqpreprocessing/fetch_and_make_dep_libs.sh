#!/bin/bash

HTSLIBVERSION="1.16"

wget https://github.com/khajoue2/libStatGen/archive/refs/tags/v1.0.15.broad.tar.gz && \
tar -zxvf v1.0.15.broad.tar.gz && \
mv libStatGen-1.0.15.broad libStatGen && \
make -C libStatGen && \

wget "https://github.com/samtools/htslib/releases/download/$HTSLIBVERSION/htslib-$HTSLIBVERSION.tar.bz2" && \
tar -jxvf "htslib-$HTSLIBVERSION.tar.bz2" && \
mv "htslib-$HTSLIBVERSION" htslib && \
cd htslib && \
./configure --disable-libcurl && \
make && \
cd .. && \

wget "https://github.com/google/googletest/archive/refs/tags/v1.13.0.tar.gz" && \
tar xf v1.13.0.tar.gz && \
mkdir -p gtest/include gtest/src && \
mv googletest-1.13.0/googletest/include/gtest gtest/include/ && \
mv googletest-1.13.0/googlemock/include/gmock gtest/include/ && \
mv googletest-1.13.0/googletest/src/* gtest/src/ && \
mv googletest-1.13.0/googlemock/src/* gtest/src/ && \
rm -f gtest/src/gmock_main.cc && \

wget http://www.cs.unc.edu/Research/compgeom/gzstream/gzstream.tgz && \
tar -xvf gzstream.tgz && \
make -C gzstream && \
echo "" && \
echo "" && \
echo "" && \
echo "If you are reading this, then the preparations succeeded." && \
echo "You should now be able to run make in this directory."

mkdir -p bin obj
