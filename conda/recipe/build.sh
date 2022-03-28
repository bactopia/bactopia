#!/bin/bash
BACTOPIA_VERSION="${PKG_VERSION%.*}.x"
BACTOPIA="${PREFIX}/share/${PKG_NAME}-${BACTOPIA_VERSION}"
mkdir -p ${PREFIX}/bin ${BACTOPIA}

chmod 775 bin/*.py
cp bin/*.py ${PREFIX}/bin

chmod 775 bin/bactopia/*
cp bin/bactopia/* ${PREFIX}/bin

# Install bactopia-py
pip3 install .

# Move bactopia nextflow
mv bin/ conda/ conf/ data/ lib/ modules/ subworkflows/ tests/ workflows/ main.nf nextflow.config ${BACTOPIA}

