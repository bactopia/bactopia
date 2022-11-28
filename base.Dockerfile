FROM condaforge/mambaforge:22.9.0-2

LABEL base.image="condaforge/mambaforge:22.9.0-2"
LABEL software="Bactopia - Base Image"
LABEL software.version="2.2.0"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robbie.petit@gmail.com"

# Dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
  curl \
  locales \
  uuid-runtime && \
  apt-get autoclean && \
  rm -rf /var/lib/apt/lists/* && \
  echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
  locale-gen en_US.UTF-8

# Add bactopia scripts, base env to path, set locale to UTF-8, set workdir
COPY bin/*.py /usr/local/bin/
ENV PATH /opt/conda/bin:$PATH
