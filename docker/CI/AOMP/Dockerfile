FROM ubuntu:20.04
MAINTAINER d.wortmann@fz-juelich.de

RUN DEBIAN_FRONTEND="noninteractive" apt update && DEBIAN_FRONTEND="noninteractive" apt install -y\
  python3 git libxml2-dev cmake vim libxml2-utils python3-pip wget

RUN wget https://github.com/ROCm-Developer-Tools/aomp/releases/download/rel_15.0-1/aomp_Ubuntu2004_15.0-1_amd64.deb  && dpkg -i aomp_*.deb
RUN python3 -m pip install pytest masci-tools pytest-html

