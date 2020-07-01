FROM nfcore/base

LABEL version="1.4.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia qc_reads process"

COPY conda/qc_reads.yml /
RUN conda env create -f qc_reads.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-qc_reads/bin:$PATH
