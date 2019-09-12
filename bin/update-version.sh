#! /bin/bash
# Updates the version numbers across the Bactopia project.
# If no user input, print usage

if [[ $# == 0 ]]; then
    echo ""
    echo "update-version.sh BACTOPIA_DIRECTORY OLD_VERSION NEW_VERSION"
    echo ""
    echo "Example Command"
    echo "update-version.sh /home/bactopia/bactopia 1.0.0 1.0.1"
    echo ""
    exit
fi


DIRECTORY=$1
OLD_VERSION=$2
NEW_VERSION=$3
if [ -z  ${DIRECTORY} ] || [ -z  ${OLD_VERSION} ] || [ -z  ${NEW_VERSION} ]; then
    echo "Got ${#} arguement"
    echo "Must give a directory, old version and new version"
    exit 1
fi

SED_CMD="echo sed -i"
if [ "$4" == "1" ]; then
    echo "In-Place edits ENABLED"
    SED_CMD="sed -i"
else
    echo "In-Place edits DISABLED (e.g. no changes will be made)"
fi

# Test $DIRECTORY points to bactopia repo
/bin/bash ${DIRECTORY}/bactopia 1> /dev/null 2> /dev/null
if [ $? -eq 0 ]; then
    # It is! Now update versions
    # Conda Files
    for file in $(ls ${DIRECTORY}/conda/*.yml); do
        ${SED_CMD} -r 's=version: '"${OLD_VERSION}"'$=version: '"${NEW_VERSION}"'=' ${file}
    done
    ${SED_CMD} 's/VERSION='"${OLD_VERSION}"'/VERSION='"${NEW_VERSION}"'/' ${DIRECTORY}/conda/README.md

    # Conda Recipe
    ${SED_CMD} 's/version = "'"${OLD_VERSION}"'"/version = "'"${NEW_VERSION}"'"/' ${DIRECTORY}/conda/anaconda/meta.yaml

    # Python Scripts
    for file in $(ls ${DIRECTORY}/bin/*.py); do
        ${SED_CMD} -r 's/VERSION = "'"${OLD_VERSION}"'"/VERSION = "'"${NEW_VERSION}"'"/' ${file}
    done

    # Docker/Singularity
    for file in $(ls ${DIRECTORY}/docker/*.Dockerfile); do
        ${SED_CMD} -r 's/version="'"${OLD_VERSION}"'"/version="'"${NEW_VERSION}"'"/' ${file}
    done
    ${SED_CMD} -r 's=container(.*):'"${OLD_VERSION}"'=container\1:'"${NEW_VERSION}"'=' ${DIRECTORY}/conf/docker.config
    ${SED_CMD} -r 's=container(.*):'"${OLD_VERSION}"'=container\1:'"${NEW_VERSION}"'=' ${DIRECTORY}/conf/singularity.config

    # Bactopia/Nextflow
    ${SED_CMD} 's/VERSION='"${OLD_VERSION}"'/VERSION='"${NEW_VERSION}"'/' ${DIRECTORY}/bactopia
    ${SED_CMD} -r "s/version = '${OLD_VERSION}'/version = '${NEW_VERSION}'/" ${DIRECTORY}/nextflow.config
    ${SED_CMD} "s/VERSION = '${OLD_VERSION}'/VERSION = '${NEW_VERSION}'/" ${DIRECTORY}/main.nf
else
    echo "Unable to execute '${DIRECTORY}/bactopia"
    echo "Please verify '${DIRECTORY}' points to the bactopia repo."
    exit 1
fi
