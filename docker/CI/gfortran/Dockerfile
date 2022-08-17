FROM ubuntu:20.04
MAINTAINER d.wortmann@fz-juelich.de

RUN DEBIAN_FRONTEND="noninteractive" apt update && DEBIAN_FRONTEND="noninteractive" apt install -y\
  python3 git gfortran libfftw3-dev g++ libblas-dev libblas-dev libxml2-dev cmake doxygen tetex-brev dvi2ps ghostcript libopenmpi-dev libscalapack-openmpi-dev vim gdb lcov libxml2-utils python3-pip

RUN pip3 install pytest masci-tools
RUN mkdir /opt/lib
RUN mkdir /opt/include
RUN ln -s /usr/lib/x86_64-linux-gnu/liblapack.so.3 /usr/lib/x86_64-linux-gnu/liblapack.so

COPY lib/* /opt/lib/
COPY include/* /opt/include/
