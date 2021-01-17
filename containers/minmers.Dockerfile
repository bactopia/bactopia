FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia - minmers"
LABEL software.version="1.5.6"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/conda/linux/minmers.yml"
LABEL conda.md5="86d0710baa67c26fde592f972389f8b1"

COPY conda/linux/minmers.yml /
RUN conda env create -q -f minmers.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-minmers/bin:$PATH
