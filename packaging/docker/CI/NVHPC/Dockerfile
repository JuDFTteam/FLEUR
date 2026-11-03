FROM --platform=linux/amd64 nvcr.io/nvidia/nvhpc:26.9-devel-cuda_multi-ubuntu22.04 AS base
MAINTAINER d.wortmann@fz-juelich.de

RUN DEBIAN_FRONTEND="noninteractive" apt update && DEBIAN_FRONTEND="noninteractive" apt install -y\
  python3 git libfftw3-dev g++ libblas-dev libblas-dev libxml2-dev cmake  vim  libxml2-utils python3-pip

RUN python3 -m pip install pytest masci-tools pytest-html py
RUN curl --proto '=https' --tlsv1.2 -sSfL https://bencher.dev/download/install-cli.sh | sh

# Build the FLEUR dependencies from the FLEUR sources in a separate stage
FROM base AS deps
ARG FLEUR_BRANCH=develop
ARG FLEUR_REVISION
COPY build-fleur-dependencies.sh /usr/local/bin/
RUN bash /usr/local/bin/build-fleur-dependencies.sh nvidia /opt/fleur-deps

FROM base
COPY --from=deps /opt/fleur-deps /opt/fleur-deps
