Bootstrap: docker
From: nfcore/base

%labels
    MAINTAINER Robert A. Petit III <robert.petit@emory.edu>
    DESCRIPTION Singularity image for Bactopia
    VERSION 1.3.1

%environment
    PATH=/opt/conda/envs/bactopia/bin:$PATH
    export PATH

%post
    /opt/conda/bin/conda create -n bactopia -c conda-forge -c bioconda bactopia=1.3.1
    /opt/conda/bin/conda clean -a
