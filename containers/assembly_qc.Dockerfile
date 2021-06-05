FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia - assembly_qc"
LABEL software.version="1.7.1"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/conda/linux/assembly_qc.yml"
LABEL conda.md5="c06a810815004979cd4eb777a2d30dd1"

COPY conda/linux/assembly_qc.yml /
RUN conda env create -q -f assembly_qc.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-assembly_qc/bin:$PATH
