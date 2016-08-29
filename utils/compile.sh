#!/bin/bash
set -e

source /hbb_exe/activate

set -x

yum install -y --quiet git wget ntp

wget http://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.tar.gz/download --no-check-certificate
tar -xvzf boost_1_55_0.tar.gz 
cd boost_1_55_0/
./bootstrap.sh --with-libraries=graph
./b2 install

export LDFLAGS="$LDFLAGS -pthread"

pwd=/io

cd $pwd
rm -rf build
mkdir build
cd $pwd/build

cmake .. -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++
make

exit $?