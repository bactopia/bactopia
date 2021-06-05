FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia Tools - eggnog"
LABEL software.version="1.7.1"
LABEL description="Functional annotation using orthologous groups"
LABEL website="https://bactopia.github.io/bactopia-tools/eggnog/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/tools/eggnog/environment-linux.md5"
LABEL conda.md5="8b51cc07886e5f754da5d76fb07a6755"

COPY tools/eggnog/environment-linux.yml /
RUN conda env create -q -f /environment-linux.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-eggnog/bin:$PATH
