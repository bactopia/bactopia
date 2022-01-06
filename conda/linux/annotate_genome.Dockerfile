FROM nfcore/base:2.1

LABEL base.image="nfcore/base:2.1"
LABEL software="Bactopia - annotate_genome"
LABEL software.version="2.0.1"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robbie.petit@gmail.com"
LABEL conda.env="bactopia/conda/linux/annotate_genome.yml"
LABEL conda.md5="9a5d66dffdc522656a63fa3e73018d19"

COPY conda/linux/annotate_genome.yml /
RUN conda env create -q -f annotate_genome.yml && conda clean -y -a 
ENV PATH /opt/conda/envs/bactopia-annotate_genome/bin:$PATH
COPY bin/*.py /opt/conda/envs/bactopia-annotate_genome/bin/
