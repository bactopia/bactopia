FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia Tools - phyloflash"
LABEL software.version="1.5.x"
LABEL description="16s assembly, alignment and tree"
LABEL website="https://bactopia.github.io/bactopia-tools/phyloflash/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"

COPY tools/phyloflash/environment-linux.yml /
RUN conda env create -q -f /environment-linux.yml \
    && conda create -n python3 -c conda-forge python=3 \
    && conda clean -y -a 

ENV PATH /opt/conda/envs/bactopia-phyloflash/bin:/opt/conda/envs/python3/bin:$PATH
