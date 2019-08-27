FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>

LABEL version="0.0.5"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia-AP qc_reads"

COPY conda/qc_reads.yml /
RUN conda env create -f qc_reads.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-qc_reads/bin:$PATH
