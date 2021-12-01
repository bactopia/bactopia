FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia - sequence_type"
LABEL software.version="2.0.0"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/conda/linux/sequence_type.yml"
LABEL conda.md5="daa1c535180307a57bd9b0ab0f4dbeed"

COPY conda/linux/sequence_type.yml /
RUN conda env create -q -f sequence_type.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-sequence_type/bin:$PATH
