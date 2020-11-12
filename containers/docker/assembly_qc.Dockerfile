FROM nfcore/base

LABEL version="1.5.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia assembly_qc process"

COPY conda/linux/assembly_qc.yml /
RUN conda env create -f assembly_qc.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-assembly_qc/bin:$PATH
