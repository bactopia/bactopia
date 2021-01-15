FROM nfcore/base:1.12.1

LABEL version="1.5.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia qc_reads process"

COPY conda/linux/qc_reads.yml /
RUN conda env create -f -q qc_reads.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-qc_reads/bin:$PATH
