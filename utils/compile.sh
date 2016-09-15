#!/bin/bash

pwd=/io

cd $pwd
rm -rf build
mkdir build
cd $pwd/build


cmake .. -DCMAKE_C_COMPILER=gcc-4.8 -DCMAKE_CXX_COMPILER=g++-4.8
make

exit $?
