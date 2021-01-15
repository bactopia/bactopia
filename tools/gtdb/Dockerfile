FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia Tools - gtdb"
LABEL software.version="1.5.x"
LABEL description="Identify marker genes and assign taxonomic classifications"
LABEL website="https://bactopia.github.io/bactopia-tools/gtdb/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"

COPY tools/gtdb/environment-linux.yml /
RUN conda env create -q -f /environment-linux.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-gtdb/bin:$PATH
