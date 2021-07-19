FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia Tools - kleborate"
LABEL software.version="1.7.1"
LABEL description="Characterize virulence and resistance in Klebsiella"
LABEL website="https://bactopia.github.io/bactopia-tools/kleborate/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/tools/kleborate/environment-linux.yml"
LABEL conda.md5="d41d8cd98f00b204e9800998ecf8427e"

COPY tools/kleborate/environment-linux.yml /
RUN conda env create -q -f /environment-linux.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-kleborate/bin:$PATH
