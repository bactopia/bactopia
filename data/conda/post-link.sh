#!/usr/bin/env bash

# Nextflow edge release
echo "Downloading Nextflow edge release..."
wget https://github.com/nextflow-io/nextflow/releases/download/v26.03.2-edge/nextflow-26.03.2-edge-dist -O ${PREFIX}/bin/nextflow
chmod +x ${PREFIX}/bin/nextflow
