#!/bin/bash
BACTOPIA="${PREFIX}/share/${PKG_NAME}-${PKG_VERSION}"
mkdir -p ${PREFIX}/bin ${BACTOPIA}

chmod 775 bin/*.py bin/helpers/*
cp bin/*.py ${PREFIX}/bin
cp bin/helpers/* ${PREFIX}/bin

chmod 775 bin/bactopia/*
cp bin/bactopia/* ${PREFIX}/bin

# Move bactopia nextflow
mv bin/ conf/ data/ modules/ subworkflows/ tests/ workflows/ main.nf citations.yml nextflow.config nextflow_schema.json ${BACTOPIA}

# Nextflow edge release
wget https://github.com/nextflow-io/nextflow/releases/download/v26.03.2-edge/nextflow-26.03.2-edge-dist -O ${PREFIX}/bin/nextflow
chmod +x ${PREFIX}/bin/nextflow
