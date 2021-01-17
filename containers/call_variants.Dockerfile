FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia - call_variants"
LABEL software.version="1.5.6"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/conda/linux/call_variants.yml"
LABEL conda.md5="fe4abb5c095507c7bd65d97c66beb1fe"

COPY conda/linux/call_variants.yml /
RUN conda env create -q -f call_variants.yml && conda clean -y -a
RUN apt-get update && apt-get install -y locales && apt-get clean -y && \
    echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && locale-gen en_US.UTF-8
ENV PATH /opt/conda/envs/bactopia-call_variants/bin:$PATH
