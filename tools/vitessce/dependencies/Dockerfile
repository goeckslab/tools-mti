FROM --platform=linux/amd64 python:3.12-slim-bookworm

SHELL ["/bin/bash", "-c"]

RUN pip install \
    --no-cache-dir \
    --no-compile \
    vitessce[all]==3.5.1 anndata==0.11.3 zarr==2.18.4 && \
    useradd -ms /bin/bash vitessce && \
    apt update && \
    apt install -y --no-install-recommends \
    openjdk-17-jdk \
    wget \
    unzip \
    && \
    apt clean && \
    rm -rf /var/lib/apt/lists/* && \
    wget 'https://downloads.openmicroscopy.org/bio-formats/8.1.0/artifacts/bftools.zip' -O bftools.zip && \
    unzip -n -j bftools.zip "bftools/*" -d /usr/local/bin/ && \
    rm -f bftools.zip && \
    mkdir -m 777 /vitessce

USER vitessce
WORKDIR /vitessce

CMD ['bash']