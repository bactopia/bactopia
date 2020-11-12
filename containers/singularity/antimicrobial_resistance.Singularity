Bootstrap: docker
From: nfcore/base

%labels
    MAINTAINER Robert A. Petit III <robert.petit@emory.edu>
    DESCRIPTION Singularity image containing requirements for the Bactopia antimicrobial_resistance process
    VERSION 1.5.x

%environment
    PATH=/opt/conda/envs/bactopia-antimicrobial_resistance/bin:$PATH
    export PATH

%files
    conda/linux/antimicrobial_resistance.yml /

%post
    /opt/conda/bin/conda env create -f /antimicrobial_resistance.yml
    /opt/conda/bin/conda clean -a
    /opt/conda/envs/bactopia-antimicrobial_resistance/bin/amrfinder -u
