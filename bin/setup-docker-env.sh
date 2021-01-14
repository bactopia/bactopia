#! /bin/bash
# Install conda-pack
set -e
ENV=$1
IS_TOOL=${2:-"0"}
conda install -c conda-forge conda-pack

if [[ "${1}" == "bactopia" ]]; then 
  echo "Building Bactopia environment"
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
elif [[ "${IS_TOOL}" == "1" ]]; then
  echo "Building Bactopia Tool ${0} environment"
  ENV="bactopia-${1}"
  conda env create -f /bactopia/tools/${1}/environment-linux.yml
else
  echo "Building Bactopia Process ${0} environment"
  ENV="bactopia-${1}"
  conda env create -f /bactopia/conda/linux/${1}.yml
fi


# Apply conda-pack
conda-pack -n ${ENV} -o /tmp/${ENV}.tar && \
  mkdir /conda && \
  cd /conda && tar xf /tmp/${ENV}.tar && \
  rm /tmp/${ENV}.tar

# Run conda-pack
/conda/bin/conda-unpack
