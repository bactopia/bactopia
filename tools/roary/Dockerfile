FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia Tools - roary"
LABEL software.version="1.5.6"
LABEL description="Pan-genome with optional core-genome tree"
LABEL website="https://bactopia.github.io/bactopia-tools/roary/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/tools/roary/environment-linux.md5"
LABEL conda.md5="e98a996999b3f5bf5e5686668cf26609"

COPY tools/roary/environment-linux.yml /
RUN conda env create -q -f /environment-linux.yml && conda clean -y -a 
ENV PATH /opt/conda/envs/bactopia-roary/bin:$PATH
