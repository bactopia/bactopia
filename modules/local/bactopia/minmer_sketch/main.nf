nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.cpus)
PROCESS_NAME = "minmer_sketch"

process MINMER_SKETCH {
    /*
    Create minmer sketches of the input FASTQs using Mash (k=21,31) and
    Sourmash (k=21,31,51)
    */
    tag "${sample}"
    label "minmer_sketch"

    publishDir "${params.outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${PROCESS_NAME}/*"
    publishDir "${params.outdir}/${sample}/minmers", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "*.{msh,sig}"

    input:
    tuple val(sample), val(single_end), path(fq)

    output:
    path("${sample}*.{msh,sig}")
    tuple val(sample), val(single_end), path("fastqs/${sample}*.fastq.gz"), path("${sample}.sig"),emit: MINMER_QUERY
    tuple val(sample), val(single_end), path("fastqs/${sample}*.fastq.gz"), path("${sample}-k31.msh"),emit: DOWNLOAD_REFERENCES
    path "${PROCESS_NAME}/*" optional true

    shell:
    fastq = single_end ? fq[0] : "${fq[0]} ${fq[1]}"
    '''
    LOG_DIR="!{PROCESS_NAME}"
    mkdir -p ${LOG_DIR}
    echo "# Timestamp" > ${LOG_DIR}/!{PROCESS_NAME}.versions
    date --iso-8601=seconds >> ${LOG_DIR}/!{PROCESS_NAME}.versions
    echo "# Mash Version" >> ${LOG_DIR}/!{PROCESS_NAME}.versions
    mash --version >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1

    echo "# Sourmash Version" >> ${LOG_DIR}/!{PROCESS_NAME}.versions
    sourmash --version >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1

    gzip -cd !{fastq} | mash sketch -o !{sample}-k21 -k 21 -s !{params.mash_sketch} -r -I !{sample} -
    gzip -cd !{fastq} | mash sketch -o !{sample}-k31 -k 31 -s !{params.mash_sketch} -r -I !{sample} -
    sourmash sketch dna -p k=21,k=31,k=51,abund,scaled=!{params.sourmash_scale} --merge !{sample} -o !{sample}.sig !{fastq}

    # pass the FASTQs along
    mkdir -p fastqs
    if [[ -L "!{fq[0]}" ]]; then
        if [ "!{single_end}" == "false" ]; then
            # Paired-End Reads
            ln -s `readlink !{fq[0]}` fastqs/!{sample}_R1.fastq.gz
            ln -s `readlink !{fq[1]}` fastqs/!{sample}_R2.fastq.gz
        else
            # Single-End Reads
            ln -s `readlink !{fq[0]}` fastqs/!{sample}.fastq.gz
        fi
    else
        if [ "!{single_end}" == "false" ]; then
            # Paired-End Reads
            cp !{fq[0]} fastqs/!{sample}_R1.fastq.gz
            cp !{fq[1]} fastqs/!{sample}_R2.fastq.gz
        else
            # Single-End Reads
            cp  !{fq[0]} fastqs/!{sample}.fastq.gz
        fi
    fi

    if [ "!{params.skip_logs}" == "false" ]; then 
        cp .command.err ${LOG_DIR}/!{PROCESS_NAME}.err
        cp .command.out ${LOG_DIR}/!{PROCESS_NAME}.out
        cp .command.sh ${LOG_DIR}/!{PROCESS_NAME}.sh || :
        cp .command.trace ${LOG_DIR}/!{PROCESS_NAME}.trace || :
    else
        rm -rf ${LOG_DIR}/
    fi
    '''

    stub:
    """
    mkdir fastqs
    mkdir ${PROCESS_NAME}
    touch fastqs/${sample}.fastq.gz
    touch ${PROCESS_NAME}/${sample}
    touch ${sample}.sig
    touch ${sample}-k31.msh
    """
}
