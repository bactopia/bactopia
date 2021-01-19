FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia"
LABEL software.version="1.5.6"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"

COPY . /bactopia
COPY bin/check-assembly-accession.py /
RUN bash /bactopia/bin/gh-actions/setup-bactopia-env.sh && \
    /opt/conda/envs/bactopia/bin/python3 check-assembly-accession.py GCF_003431365 && \
    mv /root/.config /.config && \
    chmod -R 775 /.config && \
    ln -s /.config /root/.config

ENV PATH /opt/conda/envs/bactopia/bin:$PATH
