FROM centos:6


RUN rpm --import http://ftp.scientificlinux.org/linux/scientific/5x/x86_64/RPM-GPG-KEYs/RPM-GPG-KEY-cern

RUN yum install wget -y

RUN wget http://people.centos.org/tru/devtools-2/devtools-2.repo -O /etc/yum.repos.d/devtools-2.repo
RUN yum install devtoolset-2-gcc devtoolset-2-binutils -y
RUN yum install devtoolset-2-gcc-c++ devtoolset-2-gcc-gfortran -y
RUN source /opt/rh/devtoolset-2/enable

ENV PATH=$PATH:/opt/rh/devtoolset-2/root/usr/bin/

RUN wget http://www.cmake.org/files/v3.2/cmake-3.2.2.tar.gz --no-check-certificate && tar xf cmake-3.2.2.tar.gz
RUN cd cmake-3.2.2 && ./configure && make && make install

RUN wget http://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.tar.gz --no-check-certificate
RUN tar -xvzf boost_1_55_0.tar.gz
WORKDIR /boost_1_55_0/
RUN ./bootstrap.sh --with-libraries=graph
RUN ./b2 install

RUN yum install zlib-devel -y

RUN ln -s /opt/rh/devtoolset-2/root/usr/bin/gcc /usr/bin/gcc-4.8
RUN ln -s /opt/rh/devtoolset-2/root/usr/bin/g++ /usr/bin/g++-4.8
ADD . /hinge/
WORKDIR /hinge/
RUN ./utils/build.sh
