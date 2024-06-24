FROM intel/oneapi-hpckit
MAINTAINER d.wortmann@fz-juelich.de

RUN apt-get update && apt-get -y install libxml2-dev xxd cmake python3-pip
RUN  /usr/bin/python3 -m pip install pytest masci-tools pytest-html py

RUN curl --proto '=https' --tlsv1.2 -sSfL https://bencher.dev/download/install-cli.sh | sh
#install maqao
ADD https://www.maqao.org/maqao_archive/maqao.x86_64.2.19.0.tar.xz /tmp/maqao.x86_64.2.19.0.tar.xz
RUN cd /tmp; tar xJf maqao.x86_64.2.19.0.tar.xz ; mv /tmp/maqao.x86_64.2.19.0/bin/maqao /usr/bin; rm -rf /tmp/ma*

#RUN alternatives --set python /usr/bin/python3
#glibc-static libxml2-static zlib-static vim-common docker
#COPY lib64/* /lib64/
#COPY include/* /usr/include/
ENV PATH="$PATH:/root/.cargo/bin/"