FROM nfcore/base

LABEL version="1.2.4"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia antimicrobial_resistance process"

COPY conda/antimicrobial_resistance.yml /
RUN conda env create -f antimicrobial_resistance.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-antimicrobial_resistance/bin:$PATH

RUN amrfinder -u
