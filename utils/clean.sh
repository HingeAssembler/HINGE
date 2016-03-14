#!/bin/bash

pwd=$PWD
cd $pwd/thirdparty/DAZZ_DB
make clean

cd $pwd/thirdparty//DALIGNER
make clean

cd $pwd/thirdparty/DASCRUBBER
make clean

cd $pwd/thirdparty/DEXTRACTOR
make clean

cd $pwd
mkdir build
cd $pwd/build
cmake .. -DCMAKE_C_COMPILER=gcc-4.9 -DCMAKE_CXX_COMPILER=g++-4.9
make clean

exit $?