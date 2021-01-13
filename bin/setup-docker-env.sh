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
echo "Staging Bactopia files"
BACTOPIA=/opt/conda/envs/bactopia
chmod 755 /bactopia/bactopia /bactopia/bin/helpers/*
cp /bactopia/bactopia /bactopia/bin/helpers/* ${BACTOPIA}/bin
VERSION=`${BACTOPIA}/bin/bactopia version | cut -d " " -f 2`
BACTOPIA_VERSION="${VERSION%.*}.x"
mkdir ${BACTOPIA}/share/bactopia-${BACTOPIA_VERSION}/
cp -r /bactopia/bin/ /bactopia/conda/ /bactopia/conf/ /bactopia/docs/ /bactopia/templates/ \
      /bactopia/tools/ /bactopia/main.nf /bactopia/nextflow.config ${BACTOPIA}/share/bactopia-${BACTOPIA_VERSION}/
echo "Staging complete"
# Clean up
# conda clean --all -y
# rm -rf /bactopia

# Apply conda-pack
conda install -c conda-forge conda-pack
conda-pack -n bactopia -o /tmp/bactopia.tar && \
  mkdir /conda && cd /conda && tar xf /tmp/bactopia.tar && \
  rm /tmp/bactopia.tar
  
/conda/bin/conda-unpack
