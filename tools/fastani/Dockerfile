FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia Tools - fastani"
LABEL software.version="1.7.1"
LABEL description="Pairwise average nucleotide identity"
LABEL website="https://bactopia.github.io/bactopia-tools/fastani/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/tools/fastani/environment-linux.md5"
LABEL conda.md5="d3504728b8661ade30a4b8458021b1a7"

COPY tools/fastani/environment-linux.yml /
RUN conda env create -q -f /environment-linux.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-fastani/bin:$PATH
