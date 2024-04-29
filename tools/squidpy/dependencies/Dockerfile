FROM ubuntu:jammy

ARG squidpy_version='1.4.1'
ARG python_version='3.10.*'

RUN apt update && \
    apt install -y \
        --no-install-recommends \
        python3=$python_version \
        python3-pip && \
    rm -rf /var/lib/apt/lists/* && \
    ln -s /usr/bin/python3 /usr/bin/python

RUN pip install \
        --no-cache-dir \
        --no-compile \
        squidpy==$squidpy_version

ENV \
    TMPDIR='/tmp/' \
    MPLCONFIGDIR='/tmp/mpl/' \
    NUMBA_CACHE_DIR="/tmp/numba/"