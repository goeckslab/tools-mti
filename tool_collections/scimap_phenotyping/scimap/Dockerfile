FROM ubuntu:20.04

ENV VERSION=0.17.7
RUN apt-get update && \
    apt-get install -y git python3 python3-pip unzip wget && \
    ln -s /usr/bin/python3 /usr/bin/python

RUN pip install scimap==$VERSION
