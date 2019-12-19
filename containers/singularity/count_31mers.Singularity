Bootstrap: docker
From: nfcore/base

%labels
    MAINTAINER Robert A. Petit III <robert.petit@emory.edu>
    DESCRIPTION Singularity image containing requirements for the Bactopia count_31mers process
    VERSION 1.2.3

%environment
    PATH=/opt/conda/envs/bactopia-count_31mers/bin:$PATH
    export PATH

%files
    conda/count_31mers.yml /

%post
    /opt/conda/bin/conda env create -f /count_31mers.yml
    /opt/conda/bin/conda clean -a
