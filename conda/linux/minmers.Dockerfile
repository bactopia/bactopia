FROM bactopia/base:2.2.0

LABEL base.image="bactopia/base:2.2.0"
LABEL software="Bactopia - minmers"
LABEL software.version="2.2.0"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/conda/linux/minmers.yml"
LABEL conda.md5="c92a6e86604ec915f3a6f28a40181aed"

# Bactopia environment
COPY conda/linux/minmers.yml /
RUN mamba env create -q -f minmers.yml && \
    mamba clean -a -y

# Add bactopia env to path
ENV PATH /opt/conda/envs/bactopia-minmers/bin:$PATH
