FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>

LABEL version="0.0.5"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia-AP annotate_genome"

COPY conda/annotate_genome.yml /
RUN conda env create -f annotate_genome.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-annotate_genome/bin:$PATH

RUN cd /tmp/ && \
    curl -s -o linux64.tbl2asn.gz ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz && \
    gunzip linux64.tbl2asn.gz && \
    chmod 775 linux64.tbl2asn && \
    mv linux64.tbl2asn /opt/conda/envs/bactopia-annotate_genome/bin/tbl2asn
