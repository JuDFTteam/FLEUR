FROM ubuntu
MAINTAINER d.wortmann@fz-juelich.de

RUN DEBIAN_FRONTEND="noninteractive" apt update && DEBIAN_FRONTEND="noninteractive" apt install -y\
  python3 git libxml2-dev cmake vim libxml2-utils python3-pip wget

RUN wget https://github.com/ROCm-Developer-Tools/aomp/releases/download/rel_14.0-0/aomp_Ubuntu2004_14.0-0_amd64.deb && dpkg -i aomp_Ubuntu2004_14.0-0_amd64.deb
RUN python3 -m pip install pytest 'masci-tools==0.4.8-dev5'

