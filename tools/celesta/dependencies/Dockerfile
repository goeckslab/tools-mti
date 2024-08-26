FROM jupyter/base-notebook:latest

RUN conda config --add channels bioconda \
    && mamba install -y \
    pandas \
    anndata \
    r-optparse \
    r-janitor \
    bioconductor-rhdf5 \
    r-devtools \
    r-rmixmod \
    r-spdep \
    r-reshape2 \
    r-zeallot \
    r-colorspace \
    r-viridisLite \
    r-RColorBrewer \
    r-munsell \
    r-labeling \
    r-farver \
    r-nlme \
    r-scales \
    r-mgcv \
    r-isoband \
    r-gtable \
    r-sf \
    r-s2 \
    r-rcppeigen \
    r-ggplot2 \
    && mamba clean -yaf

RUN R -e "devtools::install_github('plevritis/CELESTA')" \
    && R -e "devtools::install_github('scverse/anndataR')"