FROM nfcore/base:1.12.1

LABEL version="1.5.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia antimicrobial_resistance process"

COPY conda/linux/antimicrobial_resistance.yml /
RUN conda env create -q -f antimicrobial_resistance.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-antimicrobial_resistance/bin:$PATH

RUN amrfinder -u
