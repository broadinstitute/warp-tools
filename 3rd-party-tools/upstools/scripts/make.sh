#!/bin/sh

#  make.sh
#  upstoools
#  

if g++ -std=c++11 -pthread main.cpp cxstring.cpp trimfq.cpp sepType_DPT.cpp -o upstools
then
    echo make finished
fi

