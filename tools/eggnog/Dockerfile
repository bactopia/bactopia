FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia Tools - eggnog"
LABEL software.version="1.5.x"
LABEL description="Functional annotation using orthologous groups"
LABEL website="https://bactopia.github.io/bactopia-tools/eggnog/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"

COPY tools/eggnog/environment-linux.yml /
RUN conda env create -q -f /environment-linux.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-eggnog/bin:$PATH
