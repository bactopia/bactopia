FROM nfcore/base

LABEL version="1.5.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia Tool phyloflash"

COPY tools/phyloflash/environment-linux.yml /
RUN conda env create -f /environment-linux.yml \
    && conda create -n python3 -c conda-forge python=3 \
    && conda clean -a 

ENV PATH /opt/conda/envs/bactopia-phyloflash/bin:/opt/conda/envs/python3/bin:$PATH
