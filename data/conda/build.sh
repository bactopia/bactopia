#!/bin/bash
BACTOPIA="${PREFIX}/share/${PKG_NAME}-${PKG_VERSION}"
mkdir -p ${PREFIX}/bin ${BACTOPIA}

chmod 775 bin/*
cp bin/* ${PREFIX}/bin

# Move bactopia nextflow
mv bin/ conf/ data/ modules/ subworkflows/ tests/ workflows/ main.nf catalog.json nextflow.config nextflow_schema.json ${BACTOPIA}

# Setup the Bactopia env variables
mkdir -p ${PREFIX}/etc/conda/activate.d ${PREFIX}/etc/conda/deactivate.d
echo "export NXF_VER=26.04.0" > ${PREFIX}/etc/conda/activate.d/bactopia.sh
chmod a+x ${PREFIX}/etc/conda/activate.d/bactopia.sh

# Unset them
echo "unset NXF_VER" > ${PREFIX}/etc/conda/deactivate.d/bactopia.sh
chmod a+x ${PREFIX}/etc/conda/deactivate.d/bactopia.sh
