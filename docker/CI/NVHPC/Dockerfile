FROM nvcr.io/nvidia/nvhpc:21.2-devel-cuda11.2-ubuntu20.04
MAINTAINER d.wortmann@fz-juelich.de

RUN DEBIAN_FRONTEND="noninteractive" apt update && DEBIAN_FRONTEND="noninteractive" apt install -y\
  python3 git libfftw3-dev g++ libblas-dev libblas-dev libxml2-dev cmake  vim  libxml2-utils python3-pip

RUN python3 -m pip install pytest 'masci-tools==0.4.8-dev5'

