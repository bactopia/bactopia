FROM nfcore/base:2.1

LABEL base.image="nfcore/base:2.1"
LABEL software="Bactopia - gather_samples"
LABEL software.version="2.0.0"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robbie.petit@gmail.com"
LABEL conda.env="bactopia/conda/linux/gather_samples.yml"
LABEL conda.md5="697c00ea6e99ce3b8eba22096f819642"

COPY conda/linux/gather_samples.yml /
RUN conda env create -q -f gather_samples.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-gather_samples/bin:$PATH
