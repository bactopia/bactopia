#!/usr/bin/env bash
# update-docker
#
# Automate the building of Bactopia related Docker containers
set -e
BACTOPIA_DIR=${1:-"./"}
TYPE=${2:-""}
NAME=${3:-""}
REPOSITORY=$4:-""}
PRUNE=${5:-"0"}
VERSION=1.5.6
CONTAINER_VERSION="${VERSION%.*}.x"

function docker_push_latest {
    image=$1
    latest=$2
    has_latest=${3:-0}

    if [[ "${has_latest}" != "0" ]]; then
        echo "Pushing ${latest}"
        docker tag ${image} ${latest}
        #docker push ${latest}
    fi
}

function docker_build {
    recipe=$1
    image=$2
    latest=${3:-0}

    echo "Working on ${image}"
    docker build --rm -t ${image} -f ${recipe} .

    # Push to DockerHub
    echo "Pushing ${image}"
    #docker push ${image}
    docker_push_latest ${image} ${latest} ${latest}

    # Push to optional repos
    for repo in ${REPOSITORY}; do 
        echo "Pushing ${repo}/${image}"
        docker tag ${image} ${repo}/${image}
        #docker push ${repo}/${image}
        docker_push_latest ${image} ${repo}/${latest} ${latest}
    done

    if [[ "${PRUNE}" == "1" ]]; then
        echo "Pruning Docker Cache"
        docker image prune -a -f
        df -h
    fi
}


if [[ "${TYPE}" == "bactopia" ]]; then
    # Build Bactopia Container
    docker_build ${BACTOPIA_DIR}/Dockerfile bactopia/bactopia:${VERSION} bactopia/bactopia:latest
elif [[ "${TYPE}" == "process" ]]; then
    # Build Process Containers
    recipe_path="${BACTOPIA_DIR}/containers/docker/${NAME}.Dockerfile"
    recipe_image="bactopia/${NAME}:${CONTAINER_VERSION}"
    docker_build ${recipe_path} ${recipe_image}
elif [[ "${TYPE}" == "tools" ]]; then
    # Build Bactopia Tools containers
    recipe_path="${BACTOPIA_DIR}/tools/${NAME}"
    if [ -f "${recipe_path}/environment-linux.yml" ]; then
        docker_file="${recipe_path}/Dockerfile"
        docker_image="bactopia/tools-${NAME}:${CONTAINER_VERSION}"
        docker_build ${docker_file} ${docker_image}
    fi
else
    # Build Bactopia Tools containers
    echo "Unknown type (${TYPE}), exiting..."
fi
