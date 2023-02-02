#!/bin/bash

if [ -d gtest ]; then
  echo "gtest directory already present. If `make test` doesn't work, delete gtest and run this script again."
else
  wget "https://github.com/google/googletest/archive/refs/tags/v1.13.0.tar.gz" && \
  tar xf v1.13.0.tar.gz && \
  mkdir -p gtest/include gtest/src && \
  mv googletest-1.13.0/googletest/include/gtest gtest/include/ && \
  mv googletest-1.13.0/googlemock/include/gmock gtest/include/ && \
  mv googletest-1.13.0/googletest/src/* gtest/src/ && \
  mv googletest-1.13.0/googlemock/src/* gtest/src/ && \
  rm -f gtest/src/gmock_main.cc && \
  echo "gtest acquired. Running `make test` in this directory should now work."
fi
