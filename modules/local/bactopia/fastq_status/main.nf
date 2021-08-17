nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.cpus)
PROCESS_NAME = "fastq_status"

process FASTQ_STATUS {
    /* Determine if FASTQs are PE or SE, and if they meet minimum basepair/read counts. */
    publishDir "${params.outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${PROCESS_NAME}/*"
    publishDir "${params.outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: '*.txt'

    tag "${sample}"
    label "fastq_status"

    input:
    tuple val(sample), val(sample_type), val(single_end), path(fq), path(extra)
    
    output:
    file "*-error.txt" optional true
    tuple val(sample), val(sample_type), val(single_end), 
        path("fastqs/${sample}*.fastq.gz"), path(extra),emit: ESTIMATE_GENOME_SIZE, optional: true
    file "${PROCESS_NAME}/*" optional true

    shell:
    single_end = fq[1] == null ? true : false
    qin = sample_type.startsWith('assembly') ? 'qin=33' : 'qin=auto'
    '''
    LOG_DIR="!{PROCESS_NAME}"
    ERROR=0
    mkdir -p ${LOG_DIR}
    echo "# Timestamp" > ${LOG_DIR}/!{PROCESS_NAME}.versions
    date --iso-8601=seconds >> ${LOG_DIR}/!{PROCESS_NAME}.versions

    if [ "!{params.skip_fastq_check}" == "false" ]; then
        # Not completely sure about the inputs, so make sure they meet minimum requirements
        echo "# fastq-scan Version" >> ${LOG_DIR}/!{PROCESS_NAME}.versions
        fastq-scan -v >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1

        # Check paired-end reads have same read counts
        gzip -cd !{fq[0]} | fastq-scan > r1.json
        OPTS="--sample !{sample} --min_basepairs !{params.min_basepairs} --min_reads !{params.min_reads} --min_proportion !{params.min_proportion}"
        if [ "!{single_end}" == "false" ]; then
            if ! reformat.sh in1=!{fq[0]} in2=!{fq[1]} !{qin} out=/dev/null 2> !{sample}-paired-end-error.txt; then
                ERROR=1
                echo "!{sample} FASTQs contains an error. Please check the input FASTQs.
                    Further analysis is discontinued." | \
                sed 's/^\\s*//' >> !{sample}-paired-end-error.txt
            else
                rm -f !{sample}-paired-end-error.txt
            fi
            gzip -cd !{fq[1]} | fastq-scan > r2.json

            if ! check-fastqs.py --fq1 r1.json --fq2 r2.json ${OPTS}; then
                ERROR=1
            fi
            rm r1.json r2.json
        else
            if ! check-fastqs.py --fq1 r1.json ${OPTS}; then
                ERROR=1
            fi
            rm r1.json
        fi
    fi

    if [ "${ERROR}" -eq "0" ]; then
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
    mkdir ${PROCESS_NAME}
    mkdir fastqs
    touch ${sample}-error.txt
    touch fastqs/${sample}.fastq.gz
    touch ${PROCESS_NAME}/${sample}
    """
}
