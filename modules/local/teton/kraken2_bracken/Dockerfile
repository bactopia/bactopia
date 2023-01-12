FROM bactopia/base:2.2.0

LABEL base.image="bactopia/base:2.2.0"
LABEL software="Bactopia - kraken2_bracken"
LABEL software.version="2.2.0"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robbie.petit@gmail.com"

# kraken2_bracken environment
RUN mamba create -y -n kraken2_bracken \
    -c conda-forge \
    -c bioconda \
    bracken=2.7 \
    fastq-scan=1.0.1 \
    kraken2=2.1.2 \
    pandas=1.5.2 \
    pigz=2.6 \
    python=3.10 && \
    mamba clean -a -y && \
    mamba env export -n kraken2_bracken | grep -v "^prefix:" > conda.yml

# Add scripts and kraken2_bracken env to path
COPY ./subworkflows/local/bracken/bin/* /opt/conda/envs/kraken2_bracken/bin
ENV PATH /opt/conda/envs/kraken2_bracken/bin:$PATH
