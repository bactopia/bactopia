FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia Tools - pirate"
LABEL software.version="1.5.6"
LABEL description="Pan-genome with optional core-genome tree"
LABEL website="https://bactopia.github.io/bactopia-tools/pirate/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/tools/pirate/environment-linux.md5"
LABEL conda.md5="20ff91923ee48fef9b49d336f4603874"


COPY tools/pirate/environment-linux.yml /
RUN conda env create -q -f /environment-linux.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-pirate/bin:$PATH
