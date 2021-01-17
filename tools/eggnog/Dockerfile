FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia Tools - eggnog"
LABEL software.version="1.5.6"
LABEL description="Functional annotation using orthologous groups"
LABEL website="https://bactopia.github.io/bactopia-tools/eggnog/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/tools/eggnog/environment-linux.md5"
LABEL conda.md5="0b3e5be7cd6e13290018d3faf3858d9c"

COPY tools/eggnog/environment-linux.yml /
RUN conda env create -q -f /environment-linux.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-eggnog/bin:$PATH
