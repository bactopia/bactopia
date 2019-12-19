FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>

LABEL version="1.2.3"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia qc_reads process"

COPY conda/qc_reads.yml /
RUN conda env create -f qc_reads.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-qc_reads/bin:$PATH
