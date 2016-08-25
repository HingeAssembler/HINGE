#!/bin/bash
set -e

source /hbb_exe/activate

set -x


pwd=/io
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
cmake .. -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++
make -j 8

exit $?
