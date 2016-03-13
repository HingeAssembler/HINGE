#!/bin/bash

pwd=$PWD
cd $pwd/thirdparty/DAZZ_DB
make -j 8

cd $pwd/thirdparty//DALIGNER
make -j 8

cd $pwd/thirdparty/DASCRUBBER
make -j 8

cd $pwd/thirdparty/DEXTRACTOR
make -j 8

cd $pwd
mkdir build
cd $pwd/build
cmake .. -DCMAKE_C_COMPILER=gcc-4.9 -DCMAKE_CXX_COMPILER=g++-4.9
make -j 8

exit $?