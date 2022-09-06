FROM python:3.7.12

ENV OPS_VERSION=3.4.1

ENV PALOM_VERSION=2022.8.2

RUN apt update -y && apt install -y openslide-tools python3-opencv

RUN cd /opt && \
    wget https://github.com/openslide/openslide/releases/download/v$OPS_VERSION/openslide-$OPS_VERSION.tar.gz && \
    tar -zxvf openslide-$OPS_VERSION.tar.gz  && \
    cd openslide-$OPS_VERSION

RUN /opt/openslide-$OPS_VERSION/configure  && \
    make && \
    make install

RUN pip install zarr==2.10.3 pandas palom==$PALOM_VERSION