FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia - count_31mers"
LABEL software.version="1.7.1"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/conda/linux/count_31mers.yml"
LABEL conda.md5="70fd868316b9447d8f0c3e42e509de51"

COPY conda/linux/count_31mers.yml /
RUN conda env create -q -f count_31mers.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-count_31mers/bin:$PATH
