FROM intel/oneapi-hpckit
MAINTAINER d.wortmann@fz-juelich.de

RUN apt-get update && apt-get -y install libxml2-dev xxd
RUN  /opt/intel/oneapi/intelpython/latest/bin/python3 -m pip install pytest masci-tools
#ADD http://cmake.org/files/v3.21/cmake-3.21.2-linux-x86_64.sh /
#Certificate failed for cmake.org, hence file must be dowloaded by hand
ADD cmake-3.21.2-linux-x86_64.sh /
RUN mkdir /opt/cmake; yes|sh /cmake-3.21.2-linux-x86_64.sh --prefix=/opt/cmake; ln -sf /opt/cmake/cmake-3.21.2-linux-x86_64/bin/cmake /usr/bin/cmake
#RUN alternatives --set python /usr/bin/python3
#glibc-static libxml2-static zlib-static vim-common docker
#COPY lib64/* /lib64/
#COPY include/* /usr/include/
