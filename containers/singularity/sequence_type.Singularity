Bootstrap: docker
From: nfcore/base

%labels
    MAINTAINER Robert A. Petit III <robert.petit@emory.edu>
    DESCRIPTION Singularity image containing requirements for the Bactopia sequence_type process
    VERSION 1.3.1

%environment
    PATH=/opt/conda/envs/bactopia-sequence_type/bin:$PATH
    export PATH

%files
    conda/sequence_type.yml /

%post
    /opt/conda/bin/conda env create -f /sequence_type.yml
    /opt/conda/bin/conda clean -a
