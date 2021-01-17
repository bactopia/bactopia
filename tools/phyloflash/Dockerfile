FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia Tools - phyloflash"
LABEL software.version="1.5.6"
LABEL description="16s assembly, alignment and tree"
LABEL website="https://bactopia.github.io/bactopia-tools/phyloflash/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/tools/phyloflash/environment-linux.md5"
LABEL conda.md5="8f1eac994c2a5738927027b68843fa92"

COPY tools/phyloflash/environment-linux.yml /
RUN conda env create -q -f /environment-linux.yml \
    && conda create -n python3 -c conda-forge python=3 \
    && conda clean -y -a 

ENV PATH /opt/conda/envs/bactopia-phyloflash/bin:/opt/conda/envs/python3/bin:$PATH
