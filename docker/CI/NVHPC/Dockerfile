FROM nvcr.io/nvidia/nvhpc:23.11-devel-cuda_multi-ubuntu22.04
MAINTAINER d.wortmann@fz-juelich.de

RUN DEBIAN_FRONTEND="noninteractive" apt update && DEBIAN_FRONTEND="noninteractive" apt install -y\
  python3 git libfftw3-dev g++ libblas-dev libblas-dev libxml2-dev cmake  vim  libxml2-utils python3-pip

RUN python3 -m pip install pytest masci-tools pytest-html

