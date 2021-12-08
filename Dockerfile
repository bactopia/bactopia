FROM nfcore/base:2.1

LABEL base.image="FROM nfcore/base:2.1"
LABEL software="Bactopia"
LABEL software.version="2.0.0"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robbie.petit@gmail.com"

COPY . /bactopia
RUN bash /bactopia/bin/gh-actions/setup-bactopia-env.sh
ENV PATH /opt/conda/envs/bactopia/bin:$PATH
