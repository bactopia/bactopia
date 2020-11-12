FROM nfcore/base

LABEL version="1.5.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia annotate_genome process"

COPY conda/linux/annotate_genome.yml /
RUN conda env create -f annotate_genome.yml && conda clean -a 

ENV PATH /opt/conda/envs/bactopia-annotate_genome/bin:$PATH
