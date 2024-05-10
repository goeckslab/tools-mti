FROM ubuntu:jammy


ARG python_version='3.10.*'
ARG vitessce_version='1.0.4'
ARG anndata_version='0.8.*'
ARG numpy_version='1.23.0'
ARG bftools_url='https://downloads.openmicroscopy.org/bio-formats/6.6.1/artifacts/bftools.zip'

RUN apt update && \
    ln -snf /usr/share/zoneinfo/$CONTAINER_TIMEZONE /etc/localtime && echo $CONTAINER_TIMEZONE > /etc/timezone && \
    apt install -y \
        --no-install-recommends \
        curl \
        git \
        openjdk-8-jdk \
        unzip \
        wget \
        python3=$python_version \
        python3-pip && \
    rm -rf /var/lib/apt/lists/* && \
    ln -s /usr/bin/python3 /usr/bin/python

RUN wget $bftools_url -O bftools.zip && \
    unzip -j bftools.zip "bftools/*" -d /usr/local/bin/ && \
    rm -f bftools.zip

RUN pip install \
    --no-cache-dir \
    --no-compile \
    vitessce==$vitessce_version \
    anndata==$anndata_version \
    numpy==$numpy_version \
    scanpy \