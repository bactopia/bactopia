FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia Tools - hicap"
LABEL software.version="1.6.1"
LABEL description=""
LABEL website="https://bactopia.github.io/bactopia-tools/hicap/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/tools/hicap/environment-linux.yml"
LABEL conda.md5="00000000000000000000000000000000"

COPY tools/hicap/environment-linux.yml /
RUN conda env create -q -f /environment-linux.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-hicap/bin:$PATH
