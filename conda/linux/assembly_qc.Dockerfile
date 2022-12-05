FROM bactopia/base:2.2.0

LABEL base.image="bactopia/base:2.2.0"
LABEL software="Bactopia - assembly_qc"
LABEL software.version="2.2.0"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robbie.petit@gmail.com"
LABEL conda.env="bactopia/conda/linux/assembly_qc.yml"
LABEL conda.md5="16cfdbc493a57ce4e6a0bcd48a8b6a2e"

# Bactopia environment
COPY conda/linux/assembly_qc.yml /
RUN mamba env create -q -f assembly_qc.yml && \
    mamba clean -a -y

# Add bactopia env to path
ENV PATH /opt/conda/envs/bactopia-assembly_qc/bin:$PATH
