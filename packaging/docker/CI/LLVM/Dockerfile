FROM --platform=linux/amd64 ubuntu:24.04
MAINTAINER d.wortmann@fz-juelich.de

RUN DEBIAN_FRONTEND="noninteractive" apt update && DEBIAN_FRONTEND="noninteractive" apt install -y\
  python3 python3-pip git cmake vim wget lsb-release gnupg software-properties-common\
  libxml2-dev libxml2-utils libblas-dev liblapack-dev libfftw3-dev

RUN wget https://apt.llvm.org/llvm.sh && DEBIAN_FRONTEND="noninteractive" bash llvm.sh 21 all && rm llvm.sh
RUN ln -s /usr/bin/flang-21 /usr/local/bin/flang && ln -s /usr/bin/clang-21 /usr/local/bin/clang && ln -s /usr/bin/clang++-21 /usr/local/bin/clang++

RUN python3 -m pip install --break-system-packages pytest masci-tools pytest-html py
