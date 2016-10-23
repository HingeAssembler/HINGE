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

cd $pwd/thirdparty/racon
make modules && make tools && make -j 8

cd $pwd
mkdir build
cd $pwd/build
cmake .. -DCMAKE_INSTALL_PREFIX=../inst -DCMAKE_C_COMPILER=gcc-4.8 -DCMAKE_CXX_COMPILER=g++-4.8
make -j 8
make install

exit $?
