FROM nfcore/base:2.1

LABEL base.image="nfcore/base:2.1"
LABEL software="Bactopia - assembly_qc"
LABEL software.version="2.1.1"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robbie.petit@gmail.com"
LABEL conda.env="bactopia/conda/linux/assembly_qc.yml"
LABEL conda.md5="e7e912c75ca6345e6b81a99fb33f44cd"

COPY conda/linux/assembly_qc.yml /
RUN conda env create -q -f assembly_qc.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-assembly_qc/bin:$PATH
COPY bin/*.py /opt/conda/envs/bactopia-assembly_qc/bin/
