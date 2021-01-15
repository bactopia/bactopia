FROM nfcore/base:1.12.1

LABEL version="1.5.6"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia"

COPY . /bactopia
RUN bash /bactopia/bin/setup-docker-env.sh bactopia
ENV PATH /opt/conda/envs/bactopia/bin:$PATH
