FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia Tools - ismapper"
LABEL software.version="1.7.1"
LABEL description="Identify positions of insertion sites"
LABEL website="https://bactopia.github.io/bactopia-tools/ismapper/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/tools/ismapper/environment-linux.md5"
LABEL conda.md5="71be8fd2812c45e8fba812e3e4db7261"

COPY tools/ismapper/environment-linux.yml /
RUN conda env create -q -f /environment-linux.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-ismapper/bin:$PATH
