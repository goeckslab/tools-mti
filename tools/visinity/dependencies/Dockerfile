FROM mambaorg/micromamba

ARG TAG
ARG PIP_NO_COMPILE="true"

USER root

RUN <<APT-INSTALL
    set -e
    apt-get update
    apt-get install -y libgl1-mesa-glx curl unzip
    rm -rf /var/lib/apt/lists/*
APT-INSTALL

RUN <<VISINITY
    set -e
    cd /tmp/
    curl -L https://github.com/labsyspharm/visinity/archive/refs/tags/v$TAG.zip -o $TAG.zip
    unzip $TAG.zip
    rm -rf $TAG.zip
    mv visinity-$TAG /visinity
    useradd -ms /bin/bash visinity
    echo -e '[global]\nno-cache-dir = True' >> /etc/pip.conf
    chmod -R 0777 /visinity
    chown -R visinity /visinity
    micromamba install -n base -y -c conda-forge python=3.9 pip=24.0 openslide-python=1.3.1
    micromamba install -n base -y -f /visinity/requirements.yml
    micromamba clean --all --force-pkgs-dirs -y
VISINITY

RUN <<ENVMAKER
    echo -e '#!/bin/bash\nmicromamba activate base\n$*' >> /bin/vis-env
    chmod 755 /bin/vis-env
ENVMAKER

USER visinity

SHELL [ "micromamba", "run", "-n", "base", "/bin/bash", "-c" ]

RUN <<PIPPERY
    set -e
    pip install --no-dependencies lightkit pycave
PIPPERY

ENV BASH_ENV='/usr/local/bin/_activate_current_env.sh'

EXPOSE 8000
WORKDIR /visinity/
CMD ["python", "run.py"]