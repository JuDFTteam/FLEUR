FROM intel/oneapi-hpckit
MAINTAINER d.wortmann@fz-juelich.de

RUN apt-get update && apt-get -y install libxml2-dev xxd cmake python3-pip
RUN  /usr/bin/python3 -m pip install pytest masci-tools pytest-html
#RUN alternatives --set python /usr/bin/python3
#glibc-static libxml2-static zlib-static vim-common docker
#COPY lib64/* /lib64/
#COPY include/* /usr/include/
