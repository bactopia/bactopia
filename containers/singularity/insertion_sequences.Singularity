Bootstrap: docker
From: nfcore/base

%labels
    MAINTAINER Robert A. Petit III <robert.petit@emory.edu>
    DESCRIPTION Singularity image containing requirements for the Bactopia insertion_sequences process
    VERSION 1.2.4

%environment
    PATH=/opt/conda/envs/bactopia-insertion_sequences/bin:$PATH
    export PATH

%files
    conda/insertion_sequences.yml /

%post
    /opt/conda/bin/conda env create -f /insertion_sequences.yml
    /opt/conda/bin/conda clean -a
