#!/usr/bin/env bash
# update-docker
#
# Automate the building of Bactopia related Docker containers
set -e
REPOSITORY="quay.io/bactopia"
VERSION=1.5.5
CONTAINER_VERSION="${VERSION%.*}.x"
BACTOPIA_DIR=${1:-"./"}

function docker_build {
    recipe=$1
    image=${REPOSITORY}/$2
    latest=${3:-0}

    echo "Working on ${image}"
    docker build --rm -t ${image} -f ${recipe} .
    docker push ${image}

    if [[ "${latest}" != "0" ]]; then
        docker tag ${image} ${REPOSITORY}/${latest}
        docker push ${REPOSITORY}/${latest}
    fi
}

# Build Bactopia Container
docker_build Dockerfile bactopia:${VERSION} bactopia:latest

# Build Process Containers
for recipe in $(ls "${BACTOPIA_DIR}/containers/docker" | grep ".Dockerfile"); do
    recipe_path="${BACTOPIA_DIR}/containers/docker/${recipe}"
    recipe_name=$(echo ${recipe} | sed 's/.Dockerfile//')
    recipe_image="${recipe_name}:${CONTAINER_VERSION}"
    docker_build ${recipe_path} ${recipe_image}
done

# Build Bactopia Tools containers
for tool in $(ls "${BACTOPIA_DIR}/tools"); do
    recipe_path="${BACTOPIA_DIR}/tools/${tool}"
    if [ -f "${BACTOPIA_DIR}/tools/${tool}/environment-linux.yml" ]; then
        docker_file="${recipe_path}/Dockerfile"
        docker_image="tools-${tool}:${CONTAINER_VERSION}"
        docker_build ${docker_file} ${docker_image}
    fi
done
