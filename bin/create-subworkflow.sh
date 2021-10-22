#!/usr/bin/env bash
# build-containers
#
# Create a blank tool.
VERSION=1.7.1

if [[ $# == 0 ]]; then
    echo ""
    echo "create-tool.sh BACTOPIA_DIR SUBWORKFLOW_NAME SUBWORKFLOW_DESCRIPTION"
    echo ""
    echo "Example Command"
    echo "create-tool.sh /home/bactopia/bactopia roary 'Create a pan-genome with Roary and an optional core-genome phylogeny with IQTree.' "
    echo ""
    exit
fi

BACTOPIA_DIR=$1
SUBWORKFLOW=$2
DESCRIPTION=$3
if [ -z "${BACTOPIA_DIR}" ] || [ -z "${SUBWORKFLOW}" ] || [ -z "${DESCRIPTION}" ]; then
    echo "Got ${#} arguement"
    echo "Must give a path to Bactopia repository, tool name and tool description."
    exit 1
fi

if [ ! -d "${BACTOPIA_DIR}/sobworkflows/local/${SUBWORKFLOW}" ]; then
    cp -r ${BACTOPIA_DIR}/.skeleton/subworkflows ${BACTOPIA_DIR}/subworkflows/local/${TOOL}
    filenames=( "main.nf" "meta.yml" "test.nf" "test.yml" )
    for filename in "${filenames[@]}"; do
        sed -i -r 's/SUBWORKFLOW_NAME/'"${SUBWORKFLOW}"'/' ${BACTOPIA_DIR}/tools/${SUBWORKFLOW}/${filename}
        sed -i -r 's/SUBWORKFLOW_DESCRIPTION/'"${DESCRIPTION}"'/' ${BACTOPIA_DIR}/tools/${SUBWORKFLOW}/${filename}
    done
else
    echo "${SUBWORKFLOW} exists already, please verify. Not going to replace, exiting..."
    exit 1
fi
