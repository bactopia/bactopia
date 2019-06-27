FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>
LABEL authors="robert.petit@emory.edu" \
    description="Container image containing requirements for the Bactopia-AP ariba_analysis"

COPY conda/ariba_analysis.yml /
RUN conda env create -f ariba_analysis.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-ariba_analysis/bin:$PATH
