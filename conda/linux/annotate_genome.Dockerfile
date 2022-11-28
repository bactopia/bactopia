FROM bactopia/base:2.2.0

LABEL base.image="bactopia/base:2.2.0"
LABEL software="Bactopia - annotate_genome"
LABEL software.version="2.2.0"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robbie.petit@gmail.com"
LABEL conda.env="bactopia/conda/linux/annotate_genome.yml"
LABEL conda.md5="2d1885785fe2a7a209de1af06eea5115"

# Bactopia environment
COPY conda/linux/annotate_genome.yml /
RUN mamba env create -q -f annotate_genome.yml && \
    mamba clean -a -y 

# Add bactopia env to path
ENV PATH /opt/conda/envs/bactopia-annotate_genome/bin:$PATH

