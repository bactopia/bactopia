FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia - download_references"
LABEL software.version="1.5.6"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/conda/linux/download_references.yml"
LABEL conda.md5="78ff1aee7af081c77c4b3127c0027b6b"

COPY conda/linux/download_references.yml /
COPY bin/check-assembly-accession.py /
RUN conda env create -q -f download_references.yml && conda clean -y -a

ENV PATH /opt/conda/envs/bactopia-download_references/bin:$PATH
