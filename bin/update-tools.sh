#! /bin/bash
# Updates the conda environment yamls for Bactopia Tools to bump to latest software versions.

if [[ $# == 0 ]]; then
    echo ""
    echo "update-tools.sh BACTOPIA_DIRECTORY VERSION IS_MAC"
    echo ""
    echo "Example Command"
    echo "update-tools.sh /home/bactopia/bactopia 1.0.0"
    echo ""
    exit
fi
CONDA_DIR="${1}/tools"
VERSION=$2
IS_MAC=0
if [ "$3" == "1" ]; then
    echo "Creating Mac OS X yamls"
    IS_MAC=1
fi

function update_environment {
    # 1: template, 2: programs, 3: conda dir, 4: version, 5: is_mac
    echo "Working on ${1}"

    YAML="${3}/${1}/environment"
    if [ "$5" == 1 ]; then
        # Mac OS
        # Have to replace Mac versions of some programs (date, sed, etc...)
        conda create --quiet -y -n bactopia-${1} ${6} -c conda-forge -c bioconda ${2} coreutils sed
        conda env export --no-builds -n bactopia-${1} | grep -v "^prefix:" > ${YAML}-osx.yml
        md5 -r ${YAML}-osx.yml | cut -d " " -f 1 > ${YAML}-osx.md5
    else
        # Linux
        conda create --quiet -y -n bactopia-${1} ${6} -c conda-forge -c bioconda ${2} 
        conda env export --no-builds -n bactopia-${1} | grep -v "^prefix:"  > ${YAML}-linux.yml
        md5sum ${YAML}-linux.yml | cut -d " " -f 1 > ${YAML}-linux.md5
        head -n 1 ${YAML}-linux.md5 | xargs -I {} sed -i -E 's/(LABEL conda.md5=")(.*)(")/\1{}\3/' ${3}/${1}/Dockerfile
    fi
    
    conda env remove -n bactopia-${1}
}

# Bactopia environments
update_environment "eggnog" "eggnog-mapper" ${CONDA_DIR} ${VERSION} ${IS_MAC}
update_environment "fastani" "fastani ncbi-genome-download rename sed" ${CONDA_DIR} ${VERSION} ${IS_MAC}
update_environment "gtdb" "gtdbtk" ${CONDA_DIR} ${VERSION} ${IS_MAC}
update_environment "hicap" "hicap" ${CONDA_DIR} ${VERSION} ${IS_MAC}
update_environment "ismapper" "ismapper" ${CONDA_DIR} ${VERSION} ${IS_MAC}
update_environment "mashtree" "mashtree ncbi-genome-download rename" ${CONDA_DIR} ${VERSION} ${IS_MAC}
update_environment "phyloflash" "phyloflash mafft iqtree pigz" ${CONDA_DIR} ${VERSION} ${IS_MAC}
update_environment "pirate" "bioconductor-ggtree clonalframeml iqtree maskrc-svg ncbi-genome-download pigz pirate prokka r-dplyr r-ggplot2 r-gridextra r-phangorn rename snp-dists tbl2asn-forever" ${CONDA_DIR} ${VERSION} ${IS_MAC}
update_environment "roary" "clonalframeml iqtree maskrc-svg ncbi-genome-download pigz prokka r-ggplot2 rename roary snp-dists tbl2asn-forever" ${CONDA_DIR} ${VERSION} ${IS_MAC}
update_environment "staph-typer" "agrvate spatyper staphopia-sccmec snippy>=4.5.0" ${CONDA_DIR} ${VERSION} ${IS_MAC}
update_environment "summary" "executor jinja2" ${CONDA_DIR} ${VERSION} ${IS_MAC}

echo "Conda Last updated: " `date` > ${CONDA_DIR}/README.md
