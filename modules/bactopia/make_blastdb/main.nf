nextflow.enable.dsl = 2

process MAKE_BLASTDB {
    /* Create a BLAST database of the assembly using BLAST */
    tag "${sample}"

    publishDir "${outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${task.process}/*"
    publishDir "${outdir}/${sample}/blast", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "blastdb/*"

    input:
    tuple val(sample), val(single_end), path(fasta)

    output:
    path("blastdb/*")
    tuple val(sample), path("blastdb/*"), emit: BLAST_DB, optional:true
    file "${task.process}/*" optional true

    shell:
    '''
    LOG_DIR="!{task.process}"
    mkdir blastdb
    mkdir -p ${LOG_DIR}
    echo "# Timestamp" > ${LOG_DIR}/!{task.process}.versions
    date --iso-8601=seconds >> ${LOG_DIR}/!{task.process}.versions
    echo "# makeblastdb Version" >> ${LOG_DIR}/!{task.process}.versions
    makeblastdb -version >> ${LOG_DIR}/!{task.process}.versions 2>&1

    # Verify AWS files were staged
    if [[ ! -L "!{fasta}" ]]; then
        check-staging.py --assembly !{fasta}
    fi

    if [[ !{params.compress} == "true" ]]; then
        gzip -cd !{fasta} | \
        makeblastdb -dbtype "nucl" -title "Assembled contigs for !{sample}" -out blastdb/!{sample}
    else
        cat !{fasta} | \
        makeblastdb -dbtype "nucl" -title "Assembled contigs for !{sample}" -out blastdb/!{sample}
    fi

    if [ "!{params.skip_logs}" == "false" ]; then 
        cp .command.err ${LOG_DIR}/!{task.process}.err
        cp .command.out ${LOG_DIR}/!{task.process}.out
        cp .command.sh ${LOG_DIR}/!{task.process}.sh || :
        cp .command.trace ${LOG_DIR}/!{task.process}.trace || :
    else
        rm -rf ${LOG_DIR}/
    fi
    '''

    stub:
    """
    mkdir blastdb
    mkdir ${task.process}
    touch blastdb/${sample}
    touch ${task.process}/${sample}
    """
}

//###############
//Module testing
//###############

workflow test{

    TEST_PARAMS_CH = Channel.of([
        params.sample,
        params.single_end,
        path(params.fasta)
    ])

    make_blastdb(TEST_PARAMS_CH)
}
