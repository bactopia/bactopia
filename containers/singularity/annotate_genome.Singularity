Bootstrap: docker
From: nfcore/base

%labels
    MAINTAINER Robert A. Petit III <robert.petit@emory.edu>
    DESCRIPTION Singularity image containing requirements for the Bactopia annotate_genome process
    VERSION 1.2.4

%environment
    PATH=/opt/conda/envs/bactopia-annotate_genome/bin:$PATH
    export PATH

%files
    conda/annotate_genome.yml /

%post
    /opt/conda/bin/conda env create -f /annotate_genome.yml
    /opt/conda/bin/conda clean -a
    wget --quiet -O tbl2asn.gz ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz \
    gunzip tbl2asn.gz
    chmod 755 tbl2asn
    mv tbl2asn /usr/local/bin/tbl2asn /opt/conda/envs/bactopia-annotate_genome/bin
