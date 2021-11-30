#!/usr/bin/env bash
# update-docker
#
# Automate the building of Bactopia related Docker containers
set -e
BACTOPIA_DIR=${1:-"./"}
REPOSITORY=${2:-""}
PRUNE=${3:-"0"}
VERSION=1.7.1
CONTAINER_VERSION="${VERSION%.*}.x"

function docker_build {
    recipe=$1
    image=$2
    latest=${3:-0}

    echo "Working on ${image}"
    docker build --rm -t ${image} -f ${recipe} .

    # Push to DockerHub
    echo "Pushing ${image}"
    docker push ${image}

    if [[ "${latest}" != "0" ]]; then
        echo "Pushing ${latest}"
        docker tag ${image} ${latest}
        docker push ${latest}
    fi

    # Push to optional repos
    for repo in ${REPOSITORY}; do 
        echo "Pushing ${repo}/${image}"
        docker tag ${image} ${repo}/${image}
        docker push ${repo}/${image}

        if [[ "${latest}" != "0" ]]; then
            echo "Pushing ${repo}/${latest}"
            docker tag ${image} ${repo}/${latest}
            docker push ${repo}/${latest}
        fi
    done

    if [[ "${PRUNE}" == "1" ]]; then
        echo "Pruning Docker Cache"
        docker image prune -a -f
        df -h
    fi
}

# Build Bactopia Container
docker_build Dockerfile bactopia/bactopia:${VERSION} bactopia/bactopia:latest

# Build Process Containers
for recipe in $(ls "${BACTOPIA_DIR}/containers/docker" | grep ".Dockerfile"); do
    recipe_path="${BACTOPIA_DIR}/containers/docker/${recipe}"
    recipe_name=$(echo ${recipe} | sed 's/.Dockerfile//')
    recipe_image="bactopia/${recipe_name}:${CONTAINER_VERSION}"
    conda_yaml="${BACTOPIA_DIR}/conda/linux/${recipe}.md5"
    docker_build ${recipe_path} ${recipe_image}
done

# Build Bactopia Tools containers
for tool in $(ls "${BACTOPIA_DIR}/tools"); do
    recipe_path="${BACTOPIA_DIR}/tools/${tool}"
    if [ -f "${BACTOPIA_DIR}/tools/${tool}/environment-linux.yml" ]; then
        docker_file="${recipe_path}/Dockerfile"
        docker_image="bactopia/tools-${tool}:${CONTAINER_VERSION}"
        docker_build ${docker_file} ${docker_image}
    fi
done
