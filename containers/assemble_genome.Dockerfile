FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia - assemble_genome"
LABEL software.version="2.0.0"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/conda/linux/assemble_genome.yml"
LABEL conda.md5="558dd9e899716aa2bd0ba5977c0b0740"

COPY conda/linux/assemble_genome.yml /
RUN conda env create -q -f assemble_genome.yml && conda clean -y -a && \
    wget -O /opt/conda/envs/bactopia-assemble_genome/bin/dragonflye https://raw.githubusercontent.com/rpetit3/dragonflye/main/bin/dragonflye && \
    chmod 755 /opt/conda/envs/bactopia-assemble_genome/bin/dragonflye
ENV PATH /opt/conda/envs/bactopia-assemble_genome/bin:$PATH
