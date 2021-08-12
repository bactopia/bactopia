FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia Tools - summary"
LABEL software.version="1.5.x"
LABEL description="Generate a report summarizing Bactopia project"
LABEL website="https://bactopia.github.io/bactopia-tools/summary/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"
LABEL conda.env="bactopia/tools/summary/environment-linux.md5"
LABEL conda.md5="b51776154418c9449534eddc59cb0811"

COPY tools/summary/environment-linux.yml /
RUN conda env create -q -f /environment-linux.yml && conda clean -y -a
ENV PATH /opt/conda/envs/summary/bin:$PATH
