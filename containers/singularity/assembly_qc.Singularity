Bootstrap: docker
From: nfcore/base

%labels
    MAINTAINER Robert A. Petit III <robert.petit@emory.edu>
    DESCRIPTION Singularity image containing requirements for the Bactopia assembly_qc process
    VERSION 1.3.1

%environment
    PATH=/opt/conda/envs/bactopia-assembly_qc/bin:$PATH
    export PATH

%files
    conda/assembly_qc.yml /

%post
    /opt/conda/bin/conda env create -f /assembly_qc.yml
    /opt/conda/bin/conda clean -a
