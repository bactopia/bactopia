FROM nfcore/base:1.12.1

LABEL version="1.5.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia Too roary"

COPY tools/roary/environment-linux.yml /
RUN conda env create -q -f /environment-linux.yml && conda clean -y -a 
ENV PATH /opt/conda/envs/bactopia-roary/bin:$PATH
