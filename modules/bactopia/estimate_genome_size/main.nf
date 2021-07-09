nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.cpus)
PROCESS_NAME = "estimate_genome_size"

process ESTIMATE_GENOME_SIZE {
    /* Estimate the input genome size if not given. */
    tag "${sample}"
    label "estimate_genome_size"

    publishDir "${params.outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${PROCESS_NAME}/*"
    publishDir "${params.outdir}/${sample}", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: '*.txt'

    input:
    tuple val(sample), val(sample_type), val(single_end), path(fq), path(extra)
    val genome_size

    output:
    path "${sample}-genome-size-error.txt" optional true
    path("${sample}-genome-size.txt") optional true
    tuple val(sample), val(sample_type), val(single_end), 
        path("fastqs/${sample}*.fastq.gz"), path(extra), path("${sample}-genome-size.txt"),emit: QUALITY_CONTROL, optional: true
    path "${PROCESS_NAME}/*" optional true

    shell:
    '''
    OUTPUT="!{sample}-genome-size.txt"
    LOG_DIR="!{PROCESS_NAME}"
    mkdir -p ${LOG_DIR}
    echo "# Timestamp" > ${LOG_DIR}/!{PROCESS_NAME}.versions
    date --iso-8601=seconds >> ${LOG_DIR}/!{PROCESS_NAME}.versions

    # Verify AWS files were staged
    if [[ ! -L "!{fq[0]}" ]]; then
        if [ "!{single_end}" == "true" ]; then
            check-staging.py --fq1 !{fq[0]} --extra !{extra} --is_single
        else
            check-staging.py --fq1 !{fq[0]} --fq2 !{fq[1]} --extra !{extra}
        fi
    fi

    if [ "!{genome_size}" == "0" ]; then
        # Use mash to estimate the genome size, if a genome size cannot be
        # estimated set the genome size to 0
        echo "# Mash Version" >> ${LOG_DIR}/!{PROCESS_NAME}.versions
        mash --version >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1
        if [ "!{single_end}" == "false" ]; then
            mash sketch -o test -k 31 -m 3 -r !{fq[0]} !{fq[1]} 2>&1 | \
                grep "Estimated genome size:" | \
                awk '{if($4){printf("%d\\n", $4)}} END {if (!NR) print "0"}' > ${OUTPUT}
        else
            mash sketch -o test -k 31 -m 3 !{fq[0]} 2>&1 | \
                grep "Estimated genome size:" | \
                awk '{if($4){printf("%d\\n", $4)}} END {if (!NR) print "0"}' > ${OUTPUT}
        fi
        rm -rf test.msh
        ESTIMATED_GENOME_SIZE=`head -n1 ${OUTPUT}`

        if [ ${ESTIMATED_GENOME_SIZE} -gt "!{params.max_genome_size}" ]; then
            # Probably high coverage, try increasing number of kmer copies to 10
            if [ "!{single_end}" == "false" ]; then
                mash sketch -o test -k 31 -m 10 -r !{fq[0]} !{fq[1]} 2>&1 | \
                    grep "Estimated genome size:" | \
                    awk '{if($4){printf("%d\\n", $4)}} END {if (!NR) print "0"}' > ${OUTPUT}
            else
                mash sketch -o test -k 31 -m 10 !{fq[0]} 2>&1 | \
                    grep "Estimated genome size:" | \
                    awk '{if($4){printf("%d\\n", $4)}} END {if (!NR) print "0"}' > ${OUTPUT}
            fi
            rm -rf test.msh
        elif [ ${ESTIMATED_GENOME_SIZE} -lt "!{params.min_genome_size}" ]; then
            # Probably low coverage, try decreasing the number of kmer copies to 1
            if [ "!{single_end}" == "false" ]; then
                mash sketch -o test -k 31 -m 1 -r !{fq[0]} !{fq[1]} 2>&1 | \
                    grep "Estimated genome size:" | \
                    awk '{if($4){printf("%d\\n", $4)}} END {if (!NR) print "0"}' > ${OUTPUT}
            else
                mash sketch -o test -k 31 -m 1 !{fq[0]} 2>&1 | \
                    grep "Estimated genome size:" | \
                    awk '{if($4){printf("%d\\n", $4)}} END {if (!NR) print "0"}' > ${OUTPUT}
            fi
            rm -rf test.msh
        fi

        ESTIMATED_GENOME_SIZE=`head -n1 ${OUTPUT}`
        if [ ${ESTIMATED_GENOME_SIZE} -gt "!{params.max_genome_size}" ]; then
            rm ${OUTPUT}
            echo "!{sample} estimated genome size (${ESTIMATED_GENOME_SIZE} bp) exceeds the maximum
                    allowed genome size (!{params.max_genome_size} bp). If this is unexpected, please
                    investigate !{sample} to determine a cause (e.g. metagenomic, contaminants, etc...).
                    Otherwise, adjust the --max_genome_size parameter to fit your need. Further analysis
                    of !{sample} will be discontinued." | \
            sed 's/^\\s*//' > !{sample}-genome-size-error.txt
        elif [ ${ESTIMATED_GENOME_SIZE} -lt "!{params.min_genome_size}" ]; then
            rm ${OUTPUT}
            echo "!{sample} estimated genome size (${ESTIMATED_GENOME_SIZE} bp) is less than the minimum
                    allowed genome size (!{params.min_genome_size} bp). If this is unexpected, please
                    investigate !{sample} to determine a cause (e.g. metagenomic, contaminants, etc...).
                    Otherwise, adjust the --min_genome_size parameter to fit your need. Further analysis
                    of !{sample} will be discontinued." | \
            sed 's/^\\s*//' > !{sample}-genome-size-error.txt
        fi
    else
        # Use the genome size given by the user. (Should be >= 0)
        echo "!{genome_size}" > ${OUTPUT}
    fi

    # pass along FASTQs
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
    touch ${sample}-genome-size-error.txt
    touch ${sample}-genome-size.txt
    touch fastqs/${sample}.fastq.gz
    touch ${PROCESS_NAME}/*
    """
}
