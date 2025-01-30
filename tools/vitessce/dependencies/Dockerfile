FROM --platform=linux/amd64 python:3.12-slim-bookworm

SHELL ["/bin/bash", "-c"]

RUN pip install \
    --no-cache-dir \
    --no-compile \
    vitessce[all]==3.5.1 && \
    useradd -ms /bin/bash vitessce && \
    apt update && \
    apt install -y --no-install-recommends \
    openjdk-17-jdk \
    wget \
    unzip \
    && \
    apt clean && \
    rm -rf /var/lib/apt/lists/* && \
    wget 'https://downloads.openmicroscopy.org/bio-formats/6.6.1/artifacts/bftools.zip' -O bftools.zip && \
    unzip -n -j bftools.zip "bftools/*" -d /usr/local/bin/ && \
    rm -f bftools.zip

USER vitessce
WORKDIR /home/vitessce

CMD ['bash']