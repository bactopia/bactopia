FROM nfcore/base:2.1

LABEL base.image="nfcore/base:2.1"
LABEL software="Bactopia - sequence_type"
LABEL software.version="2.0.1"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robbie.petit@gmail.com"
LABEL conda.env="bactopia/conda/linux/sequence_type.yml"
LABEL conda.md5="db723a8a5bc44e5524db8a7993a47302"

COPY conda/linux/sequence_type.yml /
RUN conda env create -q -f sequence_type.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-sequence_type/bin:$PATH
COPY bin/*.py /opt/conda/envs/bactopia-sequence_type/bin/
