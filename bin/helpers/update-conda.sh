#! /bin/bash
# Updates the conda environment yamls to bump to latest software versions.
set -e
if [[ $# == 0 ]]; then
    echo ""
    echo "update-conda.sh BACTOPIA_DIRECTORY IS_MAC"
    echo ""
    echo "Example Command"
    echo "update-conda.sh /home/bactopia/bactopia"
    echo ""
    exit
fi


CONDA_DIR=$1/conda
IS_MAC=0
if [ "$2" == "1" ]; then
    echo "Creating Mac OS X yamls"
    CONDA_DIR="${CONDA_DIR}/mac"
    IS_MAC=1
else
    echo "Creating Linux yamls"
    CONDA_DIR="${CONDA_DIR}/linux"
fi

function update_environment {
    # 1: template, 2: programs, 3: conda dir, 4: is_mac
    echo "Working on ${1}"
    EXTRAS="coreutils sed pigz python>=3.6,<3.11"
    mamba env remove -n bactopia-${1}
    if [ "$4" == 1 ]; then
        # Mac OS
        # Have to replace Mac versions of some programs (date, sed, etc...)
        mamba create --quiet -y -n bactopia-${1} -c conda-forge -c bioconda ${2} ${EXTRAS}
        mamba env export --no-builds -n bactopia-${1} | grep -v "^prefix:" > ${3}/${1}.yml
        md5 -r ${3}/${1}.yml | cut -d " " -f 1 > ${3}/${1}.md5
    else
        # Linux
        mamba env remove -n bactopia-${1}
        mamba create --quiet -y -n bactopia-${1} -c conda-forge -c bioconda ${2} ${EXTRAS}
        mamba env export --no-builds -n bactopia-${1} | grep -v "^prefix:" > ${3}/${1}.yml
        md5sum ${3}/${1}.yml | cut -d " " -f 1 > ${3}/${1}.md5
        head -n 1 ${3}/${1}.md5 | xargs -I {} sed -i -E 's/(LABEL conda.md5=")(.*)(")/\1{}\3/' ${3}/${1}.Dockerfile
    fi
    mamba env remove -n bactopia-${1}
}

#update_environment "annotate_genome" "ncbi-amrfinderplus=3.10.45 prokka=1.14.6 tbl2asn-forever=25.7.2f" ${CONDA_DIR} ${IS_MAC}
if [ "${IS_MAC}" == "1" ]; then
#    update_environment "assemble_genome" "shovill-se=1.1.0se assembly-scan==0.4.1 unicycler=0.4.4" ${CONDA_DIR} ${IS_MAC}
#    update_environment "assembly_qc" "quast=5.2.0" ${CONDA_DIR} ${IS_MAC}
    update_environment "call_variants" "snippy=4.6.0 vcf-annotator=0.7 mash=2.3 ncbi-genome-download=0.3.1 biopython=1.77 snpeff=5.0 rename" ${CONDA_DIR} ${IS_MAC}
else
#    update_environment "assemble_genome" "shovill-se=1.1.0se dragonflye=1.0.13 nanoq=0.9.0 unicycler=0.5.0 gsl=2.6 importlib-metadata<5" ${CONDA_DIR} ${IS_MAC}
#    update_environment "assembly_qc" "checkm-genome=1.2.2 quast=5.2.0" ${CONDA_DIR} ${IS_MAC}
    update_environment "call_variants" "snippy=4.6.0 vcf-annotator=0.7 vt=2015.11.10=he941832_3 mash=2.3 ncbi-genome-download=0.3.1 biopython=1.77 snpeff=5.0 rename" ${CONDA_DIR} ${IS_MAC}
fi
#update_environment "gather_samples" "fastp=0.23.2 sra-tools=3.0.0 bbmap=39.01 art=2016.06.05 mash=2.3 ncbi-genome-download=0.3.1 fastq-dl=1.1.1 fastq-scan=1.0.1 fastqc=0.11.9 lighter=1.1.2 porechop=0.2.4 nanoq=0.9.0 nanoplot=1.40.2 rasusa=0.7.0 biopython=1.77 rename" ${CONDA_DIR} ${IS_MAC}
#update_environment "minmers" "mccortex=1.0 mash=2.3 sourmash=4.5.0" ${CONDA_DIR} ${IS_MAC}

echo "Last updated: $(date)" > ${CONDA_DIR}/README.md
