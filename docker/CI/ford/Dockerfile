FROM ubuntu:22.04
MAINTAINER d.wortmann@fz-juelich.de

RUN DEBIAN_FRONTEND="noninteractive" apt update && DEBIAN_FRONTEND="noninteractive" apt install -y\
  graphviz git python3 python3-pip python3-lxml

RUN python3 -m pip install ford

