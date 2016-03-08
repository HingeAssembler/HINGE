#!/bin/bash

pwd=$PWD
cd $pwd/DAZZ_DB
make -j 8

cd $pwd/DALIGNER
make -j 8

cd $pwd/src
mkdir build
cd $pwd/src/build
cmake .. -DCMAKE_C_COMPILER=gcc-4.9 -DCMAKE_CXX_COMPILER=g++-4.9
make -j 8

exit $?
