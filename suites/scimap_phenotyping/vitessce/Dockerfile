FROM ubuntu:20.04

RUN apt-get update && \
    ln -snf /usr/share/zoneinfo/$CONTAINER_TIMEZONE /etc/localtime && echo $CONTAINER_TIMEZONE > /etc/timezone

RUN apt-get install -y curl git openjdk-8-jdk python3 python3-pip unzip wget && \
    ln -s /usr/bin/python3 /usr/bin/python

RUN wget https://downloads.openmicroscopy.org/bio-formats/6.6.1/artifacts/bftools.zip && \
    unzip bftools.zip && \
    cp bftools/* /usr/local/bin/

RUN pip install vitessce==1.0.4 anndata argparse scanpy numpy==1.20
