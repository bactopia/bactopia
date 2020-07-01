Bootstrap: docker
From: nfcore/base

%labels
    MAINTAINER Robert A. Petit III <robert.petit@emory.edu>
    DESCRIPTION Singularity image containing requirements for the Bactopia download_references process
    VERSION 1.4.x

%environment
    PATH=/opt/conda/envs/bactopia-download_references/bin:$PATH
    export PATH

%files
    conda/download_references.yml /

%post
    /opt/conda/bin/conda env create -f /download_references.yml
    /opt/conda/bin/conda clean -a
