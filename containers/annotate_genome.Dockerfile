FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia - annotate_genome"
LABEL software.version="1.5.6"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/conda/linux/annotate_genome.yml"
LABEL conda.md5="d136f2ebe72200562c1303e49eaa8a3a"

COPY conda/linux/annotate_genome.yml /
RUN conda env create -q -f annotate_genome.yml && conda clean -y -a 
ENV PATH /opt/conda/envs/bactopia-annotate_genome/bin:$PATH
