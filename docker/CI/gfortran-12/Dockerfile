FROM  --platform=linux/amd64 ubuntu:22.04
MAINTAINER d.wortmann@fz-juelich.de

RUN DEBIAN_FRONTEND="noninteractive" apt update && DEBIAN_FRONTEND="noninteractive" apt install -y\
  python3 git gfortran-12 hdf5-tools libfftw3-dev g++ libblas-dev libblas-dev libxml2-dev cmake doxygen tetex-brev dvi2ps  libopenmpi-dev libscalapack-openmpi-dev vim gdb lcov libxml2-utils python3-pip

RUN pip3 install pytest masci-tools pytest-html
RUN mkdir /opt/lib
RUN mkdir /opt/include
RUN ln -s /usr/lib/x86_64-linux-gnu/liblapack.so.3 /usr/lib/x86_64-linux-gnu/liblapack.so
RUN update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-12 100

COPY lib/* /opt/lib/
COPY include/* /opt/include/