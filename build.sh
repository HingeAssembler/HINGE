#!/bin/bash
cd src
rm -rf build
mkdir build
cd build
cmake ..
make -j 8

