FROM ubuntu:14.04

RUN apt-get update

RUN apt-get install zlibc zlib1g zlib1g-dev -y
RUN apt-get install software-properties-common -y
RUN apt-get install build-essential wget -y
RUN apt-get install libboost-graph-dev -y 
         
RUN wget http://www.cmake.org/files/v3.2/cmake-3.2.2.tar.gz --no-check-certificate && tar xf cmake-3.2.2.tar.gz
RUN cd cmake-3.2.2 && ./configure && make && make install
         
RUN cmake --version
         
RUN add-apt-repository ppa:ubuntu-toolchain-r/test -y 
RUN apt-get update;  apt-get install gcc-4.8 g++-4.8 -y
RUN gcc-4.8 --version
RUN which gcc-4.8

ADD . /hinge/
WORKDIR /hinge/
RUN ./utils/build.sh
