Bootstrap: docker
From: nfcore/base

%labels
    MAINTAINER Robert A. Petit III <robert.petit@emory.edu>
    DESCRIPTION Singularity image containing requirements for the Bactopia ariba_analysis process
    VERSION 1.5.x

%environment
    PATH=/opt/conda/envs/bactopia-ariba_analysis/bin:$PATH
    export PATH

%files
    conda/linux/ariba_analysis.yml /

%post
    /opt/conda/bin/conda env create -f /ariba_analysis.yml
    /opt/conda/bin/conda clean -a
