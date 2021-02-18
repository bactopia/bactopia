#!/usr/bin/env bash
# build-containers
#
# Create a blank tool.
VERSION=1.6.0

if [[ $# == 0 ]]; then
    echo ""
    echo "create-tool.sh BACTOPIA_DIR TOOL_NAME TOOL_DESCRIPTION"
    echo ""
    echo "Example Command"
    echo "create-tool.sh /home/bactopia/bactopia roary 'Create a pan-genome with Roary and an optional core-genome phylogeny with IQTree.' "
    echo ""
    exit
fi

BACTOPIA_DIR=$1
TOOL=$2
DESCRIPTION=$3
if [ -z "${BACTOPIA_DIR}" ] || [ -z "${TOOL}" ] || [ -z "${DESCRIPTION}" ]; then
    echo "Got ${#} arguement"
    echo "Must give a path to Bactopia repository, tool name and tool description."
    exit 1
fi

if [ ! -d "${BACTOPIA_DIR}/tools/${TOOL}" ]; then
    cp -r ${BACTOPIA_DIR}/tools/.skeleton ${BACTOPIA_DIR}/tools/${TOOL}
    sed -i -r 's/TOOL_NAME/'"${TOOL}"'/' ${BACTOPIA_DIR}/tools/${TOOL}/Dockerfile
    sed -i -r 's/TOOL_NAME/'"${TOOL}"'/' ${BACTOPIA_DIR}/tools/${TOOL}/Singularity
    sed -i -r 's/TOOL_NAME/'"${TOOL}"'/' ${BACTOPIA_DIR}/tools/${TOOL}/nextflow.config
    sed -i -r 's/DESCRIPTION/'"${DESCRIPTION}"'/' ${BACTOPIA_DIR}/tools/${TOOL}/nextflow.config
else
    echo "${TOOL} exists already, please verify. Not going to replace, exiting..."
    exit 1
fi
