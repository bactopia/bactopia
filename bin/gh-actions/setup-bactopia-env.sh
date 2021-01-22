#! /bin/bash
# Setup Bactopia environment
# ./setup-bactopia-env.sh /path/to/conda/ /path/to/bactopia is_github_action
set -e
set -x
CONDA_PATH=${1:-"/opt/conda"}
WORK_DIR=${2:-"/bactopia"}
IS_GITHUB=${3:-"0"}
IS_GITLAB=${4:-"0"}
ENV=${5:-"bactopia"}
CONDA_CMD="create -n ${ENV}"
if [[ "${IS_GITHUB}" == "1" ]]; then
  CONDA_CMD="install"
elif [[ "${IS_GITLAB}" != "0" ]]; then
  CONDA_CMD="create --prefix ${IS_GITLAB}"
fi

# Create environment
conda ${CONDA_CMD} --quiet -y -c conda-forge -c bioconda \
  ariba \
  beautifulsoup4 \
  biopython \
  "blast>=2.10.0" \
  "bowtie2<2.4.0"  \
  cd-hit \
  conda \
  coreutils \
  executor \
  lxml \
  mash \
  ncbi-amrfinderplus \
  ncbi-genome-download \
  nextflow \
  "pysam>=0.15.3" \
  "python>3.6" \
  requests  \
  sed \
  unzip \
  wget

# Setup variables
BACTOPIA=${CONDA_PATH}/envs/${ENV}
chmod 755 ${WORK_DIR}/bactopia ${WORK_DIR}/bin/helpers/*
cp ${WORK_DIR}/bactopia ${WORK_DIR}/bin/helpers/* ${BACTOPIA}/bin
VERSION=`${BACTOPIA}/bin/bactopia version | cut -d " " -f 2`
BACTOPIA_VERSION="${VERSION%.*}.x"
BACTOPIA_SHARE="${BACTOPIA}/share/bactopia-${BACTOPIA_VERSION}/"
mkdir -p ${BACTOPIA_SHARE}

# Copy files
cp -R \
  ${WORK_DIR}/bin \
  ${WORK_DIR}/conda \
  ${WORK_DIR}/conf \
  ${WORK_DIR}/data \
  ${WORK_DIR}/templates \
  ${WORK_DIR}/tools \
  ${WORK_DIR}/main.nf \
  ${WORK_DIR}/nextflow.config \
  ${BACTOPIA_SHARE}

# Clean up
if [[ "${IS_GITHUB}" == "0" && "${IS_GITLAB}" == "0" ]]; then
  rm -rf /bactopia
  conda clean -y -a
fi
