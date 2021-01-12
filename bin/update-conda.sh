#! /bin/bash
# Updates the conda environment yamls to bump to latest software versions.

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
    # 1: template, 2: programs, 3: conda dir, 4: version, 5: is_mac, 6: extra channel
    echo "Working on ${1}"
   
    if [ "$5" == 1 ]; then
        # Mac OS
        # Have to replace Mac versions of some programs (date, sed, etc...)
        conda create -y -n bactopia-${1} ${6} -c conda-forge -c bioconda ${2} coreutils sed
        conda env export --no-builds -n bactopia-${1} > ${3}/${1}.yml
        echo "Bactopia version: ${4}" > ${3}/${1}.version
        md5 -r ${3}/${1}.yml | cut -d " " -f 1 > ${3}/${1}.md5
    else
        # Linux
        conda create -y -n bactopia-${1} ${6} -c conda-forge -c bioconda ${2} 
        conda env export --no-builds -n bactopia-${1} > ${3}/${1}.yml
        echo "Bactopia version: ${4}" > ${3}/${1}.version
        md5sum ${3}/${1}.yml | cut -d " " -f 1 > ${3}/${1}.md5
    fi
    
    conda env remove -n bactopia-${1}
}

update_environment "annotate_genome" "prokka pigz tbl2asn-forever" ${CONDA_DIR} ${VERSION} ${IS_MAC} "" 
update_environment "antimicrobial_resistance" "ncbi-amrfinderplus" ${CONDA_DIR} ${VERSION} ${IS_MAC} ""
update_environment "ariba_analysis" "ariba bowtie2=2.3.5.1" ${CONDA_DIR} ${VERSION} ${IS_MAC} ""
update_environment "assemble_genome" "shovill-se assembly-scan unicycler pigz bowtie2=2.3.5.1" ${CONDA_DIR} ${VERSION} ${IS_MAC} ""
update_environment "assembly_qc" "checkm-genome quast pigz" ${CONDA_DIR} ${VERSION} ${IS_MAC} ""
if [ "${IS_MAC}" == "1" ]; then
    update_environment "call_variants" "snippy vcf-annotator pigz vt" ${CONDA_DIR} ${VERSION} ${IS_MAC} ""
else
    update_environment "call_variants" "snippy vcf-annotator pigz vt=2015.11.10=he941832_3" ${CONDA_DIR} ${VERSION} ${IS_MAC} ""
fi
update_environment "count_31mers" "mccortex" ${CONDA_DIR} ${VERSION} ${IS_MAC} ""
update_environment "download_references" "ncbi-genome-download mash biopython python>3.6 rename" ${CONDA_DIR} ${VERSION} ${IS_MAC} ""
update_environment "gather_fastqs" "art rename ncbi-genome-download fastq-dl biopython" ${CONDA_DIR} ${VERSION} ${IS_MAC} ""
update_environment "minmers" "mash sourmash" ${CONDA_DIR} ${VERSION} ${IS_MAC} ""
update_environment "qc_reads" "bbmap fastqc fastq-scan lighter pigz" ${CONDA_DIR} ${VERSION} ${IS_MAC} ""
update_environment "sequence_type" "ariba blast bowtie2=2.3.5.1" ${CONDA_DIR} ${VERSION} ${IS_MAC} ""

echo "Last updated: " `date` > ${CONDA_DIR}/README.md
