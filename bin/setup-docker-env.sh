#! /bin/bash
# Install Bactopia environment
set -e
set -x

# Create environment
conda create --quiet -y -n bactopia -c conda-forge -c bioconda \
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

# Setup variables
BACTOPIA=/opt/conda/envs/bactopia
chmod 755 /bactopia/bactopia /bactopia/bin/helpers/*
cp /bactopia/bactopia /bactopia/bin/helpers/* ${BACTOPIA}/bin
VERSION=`${BACTOPIA}/bin/bactopia version | cut -d " " -f 2`
BACTOPIA_VERSION="${VERSION%.*}.x"
BACTOPIA_SHARE="${BACTOPIA}/share/bactopia-${BACTOPIA_VERSION}/"
mkdir -p ${BACTOPIA_SHARE}

# Copy files
cp -r \
  /bactopia/bin/ \
  /bactopia/conda/ \
  /bactopia/conf/ \
  /bactopia/docs/ \
  /bactopia/templates/ \
  /bactopia/tools/ \
  /bactopia/main.nf \
  /bactopia/nextflow.config \
  ${BACTOPIA_SHARE}

# Clean up
rm -rf /bactopia
conda clean -y -a
