FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>

LABEL version="1.2.4"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia annotate_genome process"

COPY conda/annotate_genome.yml /
RUN conda env create -f annotate_genome.yml \
    && conda clean -a \
    && wget --quiet -O tbl2asn.gz ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz \
    && gunzip tbl2asn.gz \
    && chmod 755 tbl2asn \
    && mv tbl2asn /opt/conda/envs/bactopia-annotate_genome/bin


ENV PATH /opt/conda/envs/bactopia-annotate_genome/bin:$PATH
