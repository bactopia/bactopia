FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia Tools - ectyper"
LABEL software.version="1.7.1"
LABEL description="In silico prediction of Escherichia coli serotype"
LABEL website="https://bactopia.github.io/bactopia-tools/ectyper/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/tools/ectyper/environment-linux.yml"
LABEL conda.md5="2d206941f8de75ba4d2f8beb6197db88"

COPY tools/ectyper/environment-linux.yml /
RUN conda env create -q -f /environment-linux.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-ectyper/bin:$PATH
