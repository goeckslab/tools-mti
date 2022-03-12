FROM python:3.7.12

ENV VERSION=3.4.1

RUN apt update -y && apt install -y openslide-tools python3-opencv

RUN cd /opt && \
    wget https://github.com/openslide/openslide/releases/download/v$VERSION/openslide-$VERSION.tar.gz && \
    tar -zxvf openslide-$VERSION.tar.gz  && \
    cd openslide-$VERSION

RUN /opt/openslide-$VERSION/configure  && \
    make && \
    make install

RUN pip install zarr==2.10.3 pandas palom