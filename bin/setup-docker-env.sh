#! /bin/bash

# Create the conda environment
conda create -y -n bactopia -c conda-forge -c bioconda \
  ariba \
  beautifulsoup4 \
  biopython \
  blast \
  "bowtie2<2.4.0"  \
  cd-hit \
  conda \
  coreutils \
  executor \
  lxml \
  mash \
  ncbi-genome-download \
  nextflow \
  "pysam>=0.15.3" \
  "python>3.6" \
  requests  \
  sed \
  unzip \
  wget


# Add Bactopia files
BACTOPIA=/opt/conda/envs/bactopia
cd /bactopia
chmod 755 bactopia bin/helpers/*
cp bactopia bin/helpers/* ${BACTOPIA}/bin
VERSION=`${BACTOPIA}/bin/bactopia version | cut -d " " -f 2`
BACTOPIA_VERSION="${VERSION%.*}.x"
mkdir ${BACTOPIA}/share/bactopia-${BACTOPIA_VERSION}/
cp -r bin/ conda/ conf/ docs/ templates/ tools/ main.nf nextflow.config ${BACTOPIA}/share/bactopia-${BACTOPIA_VERSION}/

# Clean up
conda clean --all -y
cd / 
rm -rf /bactopia
