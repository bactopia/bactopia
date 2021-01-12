FROM nfcore/base

LABEL version="1.5.5"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia"

COPY . /bactopia
RUN bash /bactopia/bin/setup-docker-env.sh

ENV PATH /opt/conda/envs/bactopia/bin:$PATH
