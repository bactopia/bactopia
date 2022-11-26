FROM bactopia/base:2.2.0

LABEL base.image="bactopia/base:2.2.0"
LABEL software="Bactopia - call_variants"
LABEL software.version="2.2.0"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robbie.petit@gmail.com"
LABEL conda.env="bactopia/conda/linux/call_variants.yml"
LABEL conda.md5="cfa108b3686acdd430fcc05c5c976541"

# Bactopia environment
COPY conda/linux/call_variants.yml /
RUN mamba env create -q -f call_variants.yml && \
    mamba clean -a -y

# Add bactopia env to path
ENV PATH /opt/conda/envs/bactopia-call_variants/bin:$PATH
