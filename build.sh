#!/bin/bash

pwd=$PWD
cd $pwd/DAZZ_DB
make -j 8

cd $pwd/DALIGNER
make -j 8

cd $pwd/src
mkdir build
cd $pwd/src/build
cmake ..
make -j 8

cd $pwd