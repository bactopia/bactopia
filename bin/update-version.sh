#! /bin/bash
# Updates the version numbers across the Bactopia project.
# If no user input, print usage

function generic_update {
    ${1} -r 's/'"${2}"'/'"${3}"'/' ${4}
}

function python_update {
    ${1} -r 's/VERSION = "'"${2}"'"/VERSION = "'"${3}"'"/' ${4}
}

function yaml_update {
    ${1} -r 's=version: '"${2}"'$=version: '"${3}"'=' ${4}
}

function shell_update {
    ${1} 's/VERSION='"${2}"'/VERSION='"${3}"'/' ${4}
}

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
OLD_CONTAINER="${OLD_VERSION%.*}.x"
NEW_CONTAINER="${NEW_VERSION%.*}.x"

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
    for file in $(find -not -path "*.git*" -name "*.yml" -and -not -name "mkdocs.yml" -and -not -name ".git*"); do
        yaml_update "${SED_CMD}" ${OLD_CONTAINER} ${NEW_CONTAINER} ${file}
    done

    # Python Scripts
    for file in $(find -name "*.py"); do
        python_update "${SED_CMD}" ${OLD_VERSION} ${NEW_VERSION} ${file}
    done

    # Docker/Singularity
    generic_update "${SED_CMD}" ${OLD_VERSION} ${NEW_VERSION} ./Singularity
    generic_update "${SED_CMD}" ${OLD_VERSION} ${NEW_VERSION} ./Dockerfile
    for file in $(find -mindepth 2 -not -path "*work*" -not -path "*.git*" -name "*Dockerfile" -or -name "*Singularity" -and -not -name ".git*"); do
        generic_update "${SED_CMD}" ${OLD_CONTAINER} ${NEW_CONTAINER} ${file}
    done

    # Nextflow Configs
    for file in $(find -not -path "*work*" -not -path "*.git*" -name "nextflow.config" -and -not -name ".git*"); do
        generic_update "${SED_CMD}" ${OLD_VERSION} ${NEW_VERSION} ${file}
        generic_update "${SED_CMD}" ${OLD_CONTAINER} ${NEW_CONTAINER} ${file}
    done

    # Bactopia/Nextflow
    shell_update "${SED_CMD}" ${OLD_VERSION} ${NEW_VERSION} ${DIRECTORY}/bactopia
    shell_update "${SED_CMD}" ${OLD_VERSION} ${NEW_VERSION} ${DIRECTORY}/bin/build-containers.sh
    shell_update "${SED_CMD}" ${OLD_VERSION} ${NEW_VERSION} ${DIRECTORY}/bin/create-tool.sh

else
    echo "Unable to execute '${DIRECTORY}/bactopia"
    echo "Please verify '${DIRECTORY}' points to the bactopia repo."
    exit 1
fi
