FROM ubuntu:20.04
MAINTAINER d.wortmann@fz-juelich.de

RUN DEBIAN_FRONTEND="noninteractive" apt update && DEBIAN_FRONTEND="noninteractive" apt install -y\
  python3 git gfortran libfftw3-dev g++ libblas-dev libblas-dev libxml2-dev cmake doxygen libopenmpi-dev libscalapack-openmpi-dev vim gdb lcov libxml2-utils

#RUN apt update && apt install -y \
#    gfortran-7 g++-7 gcc-7 \
#    build-essential \
#    automake libtool \
#    binutils-gold lcov \
#    gsl-bin libgsl-dev \
#    python3-pytest \
#    valgrind vim git gdb cmake libxml2-dev liblapack-dev libblas-dev doxygen libopenmpi-dev libscalapack-mpi-dev

#RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 100
#RUN update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-7 100
#RUN update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-7 100

RUN mkdir /opt/lib
RUN mkdir /opt/include
RUN ln -s /usr/lib/x86_64-linux-gnu/liblapack.so.3 /usr/lib/x86_64-linux-gnu/liblapack.so

COPY lib/* /opt/lib/
COPY include/* /opt/include/
