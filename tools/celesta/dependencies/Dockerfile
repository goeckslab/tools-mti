FROM jupyter/base-notebook:latest

COPY requirements.txt /tmp/

RUN conda config --set channel_priority strict \
    && mamba install -y --file /tmp/requirements.txt \
    && mamba clean -yaf

RUN R -e "devtools::install_github('plevritis/CELESTA')"
