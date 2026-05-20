FROM mcr.microsoft.com/devcontainers/base:ubuntu-24.04

# Core build tools and scientific dependencies for FLEUR development.
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    ninja-build \
    pkg-config \
    gcc \
    g++ \
    gfortran \
    openmpi-bin \
    libopenmpi-dev \
    libopenblas-openmp-dev \
    libscalapack-openmpi-dev \
    libfftw3-dev \
    libxml2-dev \
    libhdf5-openmpi-dev \
    python3-venv \
    python3.12-venv \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install --no-cache-dir --break-system-packages pre-commit fprettify fortls

ENV CC=gcc
ENV CXX=g++
ENV FC=gfortran
ENV MPICC=mpicc
ENV MPICXX=mpicxx
ENV MPIFC=mpifort
ENV BLA_VENDOR=OpenBLAS
