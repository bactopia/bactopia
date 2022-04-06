FROM nfcore/base:2.1

LABEL base.image="nfcore/base:2.1"
LABEL software="Bactopia - gather_samples"
LABEL software.version="2.0.3"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robbie.petit@gmail.com"
LABEL conda.env="bactopia/conda/linux/gather_samples.yml"
LABEL conda.md5="b598fdf7db22aa68db9fe001f2c42c66"

COPY conda/linux/gather_samples.yml /
RUN conda env create -q -f gather_samples.yml && conda clean -y -a
RUN apt-get update && apt-get --quiet install --yes curl uuid-runtime && apt-get clean
RUN chmod -R 777 /root && \
    mkdir -p /root/.ncbi/ && \
    mkdir /.ncbi && \
    printf '/LIBS/GUID = "%s"\n' $(uuidgen) > /root/.ncbi/user-settings.mkfg && \
    cp /root/.ncbi/user-settings.mkfg /.ncbi/user-settings.mkfg

ENV PATH /opt/conda/envs/bactopia-gather_samples/bin:$PATH
COPY bin/*.py /opt/conda/envs/bactopia-gather_samples/bin/
