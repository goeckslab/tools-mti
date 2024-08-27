FROM jupyter/base-notebook:latest

RUN conda config --add channels bioconda \
    && mamba install -y \
    pandas=1.5.3 \
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
    r-s2 \
    r-rcppeigen \
    r-ggplot2 \
    && mamba clean -yaf

RUN R -e "devtools::install_github('plevritis/CELESTA')" \
    && R -e "devtools::install_github('scverse/anndataR@220f977c15669a22207b37d437fc6d00b060b95f')" \
    && R -e "devtools::install_github('r-spatial/sf')"