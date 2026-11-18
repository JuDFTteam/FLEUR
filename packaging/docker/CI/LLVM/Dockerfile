FROM --platform=linux/amd64 ubuntu:24.04 AS base
MAINTAINER d.wortmann@fz-juelich.de

RUN DEBIAN_FRONTEND="noninteractive" apt update && DEBIAN_FRONTEND="noninteractive" apt install -y\
  python3 python3-pip python3-venv git cmake vim wget lsb-release gnupg software-properties-common\
  libxml2-dev libxml2-utils libblas-dev liblapack-dev libfftw3-dev

RUN wget https://apt.llvm.org/llvm.sh && DEBIAN_FRONTEND="noninteractive" bash llvm.sh 22 all && rm llvm.sh
RUN DEBIAN_FRONTEND="noninteractive" apt install -y flang-22
RUN ln -s /usr/bin/flang-22 /usr/local/bin/flang && ln -s /usr/bin/clang-22 /usr/local/bin/clang && ln -s /usr/bin/clang++-22 /usr/local/bin/clang++

RUN python3 -m pip install --break-system-packages pytest masci-tools pytest-html py

# Build the FLEUR dependencies from the FLEUR sources in a separate stage
FROM base AS deps
ARG FLEUR_BRANCH=develop
ARG FLEUR_REVISION
COPY build-fleur-dependencies.sh /usr/local/bin/
RUN bash /usr/local/bin/build-fleur-dependencies.sh flang /opt/fleur-deps

FROM base
COPY --from=deps /opt/fleur-deps /opt/fleur-deps
