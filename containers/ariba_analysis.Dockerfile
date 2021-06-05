FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia - ariba_analysis"
LABEL software.version="1.7.1"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/conda/linux/ariba_analysis.yml"
LABEL conda.md5="9612da6e987bbc02d7a421fabdbf6da0"

COPY conda/linux/ariba_analysis.yml /
RUN conda env create -q -f ariba_analysis.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-ariba_analysis/bin:$PATH
