FROM nfcore/base

LABEL version="1.2.4"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia ariba_analysis process"

COPY conda/ariba_analysis.yml /
RUN conda env create -f ariba_analysis.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-ariba_analysis/bin:$PATH
