#!/bin/bash

wget "https://github.com/google/googletest/archive/refs/tags/v1.13.0.tar.gz" && \
tar xf v1.13.0.tar.gz && \
mkdir -p gtest/include gtest/src && \
mv googletest-1.13.0/googletest/include/gtest gtest/include/ && \
mv googletest-1.13.0/googlemock/include/gmock gtest/include/ && \
mv googletest-1.13.0/googletest/src/* gtest/src/ && \
mv googletest-1.13.0/googlemock/src/* gtest/src/ && \
rm -f gtest/src/gmock_main.cc && \

echo "gtest acquired. Running `make test` in this directory should now work."
