Bootstrap: docker
From: nfcore/base

%labels
    MAINTAINER Robert A. Petit III <robert.petit@emory.edu>
    DESCRIPTION Singularity image containing requirements for the Bactopia annotate_genome process
    VERSION 1.2.2

%environment
    PATH=/opt/conda/envs/bactopia-annotate_genome/bin:$PATH
    export PATH

%files
    conda/annotate_genome.yml /

%post
    /opt/conda/bin/conda env create -f /annotate_genome.yml
    /opt/conda/bin/conda clean -a
