#! /bin/bash
# Updates the conda environment yamls to bump to latest software versions.
set -x
set -e
if [[ $# == 0 ]]; then
    echo ""
    echo "update-conda.sh BACTOPIA_DIRECTORY VERSION IS_MAC"
    echo ""
    echo "Example Command"
    echo "update-conda.sh /home/bactopia/bactopia 1.0.0"
    echo ""
    exit
fi


CONDA_DIR=$1/conda
DOCKER_DIR=$1/containers
VERSION=$2
IS_MAC=0
if [ "$3" == "1" ]; then
    echo "Creating Mac OS X yamls"
    CONDA_DIR="${CONDA_DIR}/mac"
    IS_MAC=1
else
    echo "Creating Linux yamls"
    CONDA_DIR="${CONDA_DIR}/linux"
fi

function update_environment {
    # 1: template, 2: programs, 3: conda dir, 4: docker dir, 5: version, 6: is_mac
    echo "Working on ${1}"
    EXTRAS="coreutils sed pigz python>=3.6"
   
    if [ "$6" == 1 ]; then
        # Mac OS
        # Have to replace Mac versions of some programs (date, sed, etc...)
        mamba create --quiet -y -n bactopia-${1} -c conda-forge -c bioconda ${2} ${EXTRAS}
        conda env export --no-builds -n bactopia-${1} | grep -v "^prefix:" > ${3}/${1}.yml
        md5 -r ${3}/${1}.yml | cut -d " " -f 1 > ${3}/${1}.md5
    else
        # Linux
        mamba create --quiet -y -n bactopia-${1} -c conda-forge -c bioconda ${2} ${EXTRAS}
        conda env export --no-builds -n bactopia-${1} | grep -v "^prefix:" > ${3}/${1}.yml
        md5sum ${3}/${1}.yml | cut -d " " -f 1 > ${3}/${1}.md5
        head -n 1 ${3}/${1}.md5 | xargs -I {} sed -i -E 's/(LABEL conda.md5=")(.*)(")/\1{}\3/' ${4}/${1}.Dockerfile
    fi
    
    conda env remove -n bactopia-${1}
}

update_environment "annotate_genome" "prokka=1.14.6 tbl2asn-forever=25.7.2f" ${CONDA_DIR} ${DOCKER_DIR} ${VERSION} ${IS_MAC}
update_environment "antimicrobial_resistance" "ncbi-amrfinderplus=3.10.1 blast=2.11.0" ${CONDA_DIR} ${DOCKER_DIR} ${VERSION} ${IS_MAC}
update_environment "ariba_analysis" "ariba=2.14.6 bowtie2=2.3.5.1 tbb=2020.2" ${CONDA_DIR} ${DOCKER_DIR} ${VERSION} ${IS_MAC}
update_environment "assemble_genome" "shovill-se=1.1.0se assembly-scan=0.3.0 unicycler=0.4.8 bowtie2=2.3.5.1" ${CONDA_DIR} ${DOCKER_DIR} ${VERSION} ${IS_MAC}
if [ "${IS_MAC}" == "1" ]; then
    update_environment "assembly_qc" "quast=5.0.2" ${CONDA_DIR} ${DOCKER_DIR} ${VERSION} ${IS_MAC}
    update_environment "call_variants" "snippy=4.6.0 vcf-annotator=0.6" ${CONDA_DIR} ${DOCKER_DIR} ${VERSION} ${IS_MAC}
else
    update_environment "assembly_qc" "checkm-genome=1.1.3 quast=5.0.2" ${CONDA_DIR} ${DOCKER_DIR} ${VERSION} ${IS_MAC}
    update_environment "call_variants" "snippy=4.6.0 vcf-annotator=0.6 vt=2015.11.10=he941832_3" ${CONDA_DIR} ${DOCKER_DIR} ${VERSION} ${IS_MAC}
fi
update_environment "count_31mers" "mccortex=1.0" ${CONDA_DIR} ${DOCKER_DIR} ${VERSION} ${IS_MAC}
update_environment "download_references" "ncbi-genome-download=0.3.0 mash=2.3 biopython rename" ${CONDA_DIR} ${DOCKER_DIR} ${VERSION} ${IS_MAC}
update_environment "gather_fastqs" "art=2016.06.05 ncbi-genome-download=0.3.0 fastq-dl=1.0.6 biopython rename" ${CONDA_DIR} ${DOCKER_DIR} ${VERSION} ${IS_MAC}
update_environment "minmers" "mash=2.3 sourmash=4.0.0" ${CONDA_DIR} ${DOCKER_DIR} ${VERSION} ${IS_MAC}
update_environment "qc_reads" "bbmap=38.90 fastqc=0.11.9 fastq-scan=0.4.3 lighter=1.1.2" ${CONDA_DIR} ${DOCKER_DIR} ${VERSION} ${IS_MAC}
update_environment "sequence_type" "ariba=2.14.6 blast=2.11.0 bowtie2=2.3.5.1 tbb=2020.2" ${CONDA_DIR} ${DOCKER_DIR} ${VERSION} ${IS_MAC}

echo "Last updated: " `date` > ${CONDA_DIR}/README.md
