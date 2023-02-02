#!/bin/bash
BACTOPIA_VERSION="${PKG_VERSION%.*}.x"
BACTOPIA="${PREFIX}/share/${PKG_NAME}-${BACTOPIA_VERSION}"
mkdir -p ${PREFIX}/bin ${BACTOPIA}

chmod 775 bin/*.py
cp bin/*.py ${PREFIX}/bin

chmod 775 bin/bactopia/*
cp bin/bactopia/* ${PREFIX}/bin

# Install bactopia-py
$PYTHON -m pip install . --no-deps --ignore-installed -vv

# Download datasets
wget --quiet -O data/datasets/amrfinderdb.tar.gz https://datasets.bactopia.com/datasets/v${PKG_VERSION}/amrfinderdb.tar.gz
wget --quiet -O data/datasets/mlst.tar.gz https://datasets.bactopia.com/datasets/v${PKG_VERSION}/mlst.tar.gz

# Move bactopia nextflow
mv bin/ conda/ conf/ data/ lib/ modules/ subworkflows/ tests/ workflows/ main.nf nextflow.config ${BACTOPIA}
