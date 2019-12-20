FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>

LABEL version="1.2.4"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia Tool TOOL_NAME"

RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/TOOL_NAME/bin:$PATH
