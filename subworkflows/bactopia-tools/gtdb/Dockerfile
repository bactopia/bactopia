FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia Tools - gtdb"
LABEL software.version="1.7.1"
LABEL description="Identify marker genes and assign taxonomic classifications"
LABEL website="https://bactopia.github.io/bactopia-tools/gtdb/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/tools/gtdb/environment-linux.md5"
LABEL conda.md5="0e84e1c6e2602f5dbb4d440e0ad27077"

COPY tools/gtdb/environment-linux.yml /
RUN conda env create -q -f /environment-linux.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-gtdb/bin:$PATH
