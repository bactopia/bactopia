FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>

LABEL version="1.2.1"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia call_variants process"

COPY conda/call_variants.yml /
RUN conda env create -f call_variants.yml && conda clean -a
RUN apt-get update && apt-get install -y locales && apt-get clean -y && \
    echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && locale-gen en_US.UTF-8
ENV PATH /opt/conda/envs/bactopia-call_variants/bin:$PATH
