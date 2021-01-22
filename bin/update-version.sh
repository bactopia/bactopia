#! /bin/bash
# Updates the version numbers across the Bactopia project.
# If no user input, print usage

function generic_update {
    ${1} -r 's/'"${2}"'/'"${3}"'/' ${4}
}

function python_update {
    ${1} -r 's/VERSION = "'"${2}"'"/VERSION = "'"${3}"'"/' ${4}
}

function conda_update {
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
    IGNORE=${DIRECTORY}/data/version-ignore.txt
    EXCLUDE=${DIRECTORY}/data/version-excludes.txt
    for file in $(find -type f | grep -v -f ${IGNORE} | xargs -I {} grep -i -H "version" {} | grep -v -f ${EXCLUDE} | cut -d ":" -f 1 | sort | uniq); do
        if [[ "${file}" == *"bactopia" ]]; then
            # bactopia
            shell_update "${SED_CMD}" ${OLD_VERSION} ${NEW_VERSION} ${file}
        elif [[ "${file}" == *".version" ]]; then
            # Conda
            conda_update "${SED_CMD}" ${OLD_CONTAINER} ${NEW_CONTAINER} ${file}
        elif [[ "${file}" == *"Dockerfile" ]]; then
            # Docker
            generic_update "${SED_CMD}" ${OLD_VERSION} ${NEW_VERSION} ${file}
        elif [[ "${file}" == *"nextflow.config" ]]; then
            # Nextflow Config
            generic_update "${SED_CMD}" ${OLD_VERSION} ${NEW_VERSION} ${file}
            generic_update "${SED_CMD}" ${OLD_CONTAINER} ${NEW_CONTAINER} ${file}
        elif [[ "${file}" == *"Singularity" ]]; then
            # Singularity
            generic_update "${SED_CMD}" ${OLD_VERSION} ${NEW_VERSION} ${file}
        elif [[ "${file}" == *".py" ]]; then
            # Python
            python_update "${SED_CMD}" ${OLD_VERSION} ${NEW_VERSION} ${file}
        elif [[ "${file}" == *".sh" ]]; then
            # Shell
            shell_update "${SED_CMD}" ${OLD_VERSION} ${NEW_VERSION} ${file}
        else
            echo "Unknown: ${file}"
        fi
    done
else
    echo "Unable to execute '${DIRECTORY}/bactopia"
    echo "Please verify '${DIRECTORY}' points to the bactopia repo."
    exit 1
fi
