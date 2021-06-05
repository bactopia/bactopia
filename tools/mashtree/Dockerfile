FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia Tools - mashtree"
LABEL software.version="1.7.1"
LABEL description="Trees based on Mash distances"
LABEL website="https://bactopia.github.io/bactopia-tools/mashtree/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/tools/mashtree/environment-linux.md5"
LABEL conda.md5="d7a5e793818190a07b28aeff0a687131"

COPY tools/mashtree/environment-linux.yml /
RUN conda env create -q -f /environment-linux.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-mashtree/bin:$PATH
