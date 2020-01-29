FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>

LABEL version="1.2.4"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia Tool: roary"

RUN conda env create -f /environment.yml && \
    && conda clean -a \
    && wget --quiet -O tbl2asn.gz ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz \
    && gunzip tbl2asn.gz \
    && chmod 755 tbl2asn \
    && mv tbl2asn /opt/conda/envs/bactopia/bin

ENV PATH /opt/conda/envs/bactopia-roary/bin:$PATH
