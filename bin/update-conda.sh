#! /bin/bash
# Updates the conda environment yamls to bump to latest software versions.
set -x

if [[ $# == 0 ]]; then
    echo ""
    echo "update-conda.sh BACTOPIA_DIRECTORY VERSION"
    echo ""
    echo "Example Command"
    echo "update-conda.sh /home/bactopia/bactopia 1.0.0"
    echo ""
    exit
fi

CONDA_DIR=$1/conda
VERSION=$2

function update_environment {
    # 1: template, 2: programs, 3: conda dir, 4: version, 5: extra channel
    echo "Working on ${1}"
    conda create -y -n bactopia-${1} ${5} -c conda-forge -c bioconda ${2} 
    conda env export -n bactopia-${1} | \
        grep -v prefix | \
        sed -r 's=channels:=version: '"${4}"'\nchannels:=' > ${3}/${1}.yml
    conda env remove -n bactopia-${1}
}

update_environment "annotate_genome" "prokka pigz tbl2asn-forever" ${CONDA_DIR} ${VERSION} ""
update_environment "antimicrobial_resistance" "ncbi-amrfinderplus" ${CONDA_DIR} ${VERSION} ""
update_environment "ariba_analysis" "ariba" ${CONDA_DIR} ${VERSION} ""
update_environment "assemble_genome" "shovill assembly-scan unicycler pigz" ${CONDA_DIR} ${VERSION} "-c rpetit3"
update_environment "assembly_qc" "checkm-genome quast pigz" ${CONDA_DIR} ${VERSION} ""
update_environment "call_variants" "snippy samtools=1.9 vcf-annotator pigz" ${CONDA_DIR} ${VERSION} ""
update_environment "count_31mers" "mccortex" ${CONDA_DIR} ${VERSION} ""
update_environment "download_references" "ncbi-genome-download mash biopython python>3.6" ${CONDA_DIR} ${VERSION} ""
update_environment "gather_fastqs" "aspera-connect art rename ncbi-genome-download fastq-dl" ${CONDA_DIR} ${VERSION} "-c rpetit3"
update_environment "minmers" "mash sourmash" ${CONDA_DIR} ${VERSION} ""
update_environment "qc_reads" "bbmap fastqc fastq-scan lighter pigz" ${CONDA_DIR} ${VERSION} ""
update_environment "sequence_type" "ariba blast" ${CONDA_DIR} ${VERSION} ""

echo "Last updated: " `date` > ${CONDA_DIR}/README.md
