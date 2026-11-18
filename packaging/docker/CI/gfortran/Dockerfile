FROM --platform=linux/amd64 ubuntu:26.10 AS base
MAINTAINER d.wortmann@fz-juelich.de

RUN DEBIAN_FRONTEND="noninteractive" apt update && DEBIAN_FRONTEND="noninteractive" apt install -y\
  curl python3 git gfortran-16 gcc-16 g++-16 hdf5-tools libhdf5-openmpi-dev libfftw3-dev libblas-dev libxml2-dev cmake tetex-brev dvi2ps  libopenmpi-dev libscalapack-openmpi-dev vim gdb lcov libxml2-utils python3-pip python3-venv valgrind

# Make GCC 16 the default compiler, also for the MPI wrappers
RUN for c in gfortran gcc g++; do ln -s /usr/bin/$c-16 /usr/local/bin/$c; done && ln -s /usr/bin/gcc-16 /usr/local/bin/cc && ln -s /usr/bin/g++-16 /usr/local/bin/c++

#RUN pip3 install --break-system-packages masci-tools
RUN ln -s /usr/lib/x86_64-linux-gnu/liblapack.so.3 /usr/lib/x86_64-linux-gnu/liblapack.so

# Build libxc, Wannier90 and ELPA from the FLEUR sources in a separate stage
FROM base AS deps
ARG FLEUR_BRANCH=develop
ARG FLEUR_REVISION
COPY build-fleur-dependencies.sh /usr/local/bin/
RUN bash /usr/local/bin/build-fleur-dependencies.sh mpi /opt/fleur-deps

FROM base
RUN curl --proto '=https' --tlsv1.2 -sSfL https://bencher.dev/download/install-cli.sh |sed -e s/musl-static/gnu/| sh
COPY --from=deps /opt/fleur-deps /opt/fleur-deps
ENV PATH="$PATH:/root/.cargo/bin/"
