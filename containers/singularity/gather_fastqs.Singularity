Bootstrap: docker
From: nfcore/base

%labels
    MAINTAINER Robert A. Petit III <robert.petit@emory.edu>
    DESCRIPTION Singularity image containing requirements for the Bactopia gather_fastqs process
    VERSION 1.5.x

%environment
    PATH=/opt/conda/envs/bactopia-gather_fastqs/bin:$PATH
    export PATH

%files
    conda/linux/gather_fastqs.yml /

%post
    /opt/conda/bin/conda env create -f /gather_fastqs.yml
    /opt/conda/bin/conda clean -a
