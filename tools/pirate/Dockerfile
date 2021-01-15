FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia Tools - pirate"
LABEL software.version="1.5.x"
LABEL description="Pan-genome with optional core-genome tree"
LABEL website="https://bactopia.github.io/bactopia-tools/pirate/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"

COPY tools/pirate/environment-linux.yml /
RUN conda env create -q -f /environment-linux.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-pirate/bin:$PATH
