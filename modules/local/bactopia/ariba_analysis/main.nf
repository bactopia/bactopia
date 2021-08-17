nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.cpus)
PROCESS_NAME = "ariba_analysis"

process ARIBA_ANALYSIS {
    /* Run reads against all available (if any) ARIBA datasets */
    tag "${sample} - ${dataset_name}"
    label "ariba_analysis"

    publishDir "${params.outdir}/${sample}/logs", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${PROCESS_NAME}/*"
    publishDir "${params.outdir}/${sample}/ariba", mode: "${params.publish_mode}", overwrite: params.overwrite, pattern: "${dataset_name}/*"

    input:
    tuple val(sample), val(single_end), path(fq)
    each path(dataset)

    output:
    file "${dataset_name}/*"
    file "${PROCESS_NAME}/*" optional true

    when:
    single_end == false

    shell:
    dataset_tarball = dataset.getName()
    dataset_name = dataset_tarball.replace('.tar.gz', '')
    spades_options = params.spades_options ? "--spades_options '${params.spades_options}'" : ""
    noclean = params.ariba_no_clean ? "--noclean" : ""
    '''
    LOG_DIR="!{PROCESS_NAME}"
    mkdir -p ${LOG_DIR}
    echo "# Timestamp" > ${LOG_DIR}/!{PROCESS_NAME}.versions
    date --iso-8601=seconds >> ${LOG_DIR}/!{PROCESS_NAME}.versions

    # Print captured STDERR incase of exit
    function print_stderr {
        cat .command.err 1>&2
        ls ${LOG_DIR}/ | grep ".err" | xargs -I {} cat ${LOG_DIR}/{} 1>&2
    }
    trap print_stderr EXIT

    tar -xzvf !{dataset_tarball}
    mv !{dataset_name} !{dataset_name}db
    # ariba Version
    echo "# Ariba Version" >> ${LOG_DIR}/!{PROCESS_NAME}.versions
    ariba version >> ${LOG_DIR}/!{PROCESS_NAME}.versions 2>&1
    ariba run !{dataset_name}db !{fq} !{dataset_name} \
            --nucmer_min_id !{params.nucmer_min_id} \
            --nucmer_min_len !{params.nucmer_min_len} \
            --nucmer_breaklen !{params.nucmer_breaklen} \
            --assembly_cov !{params.assembly_cov} \
            --min_scaff_depth !{params.min_scaff_depth} \
            --assembled_threshold !{params.assembled_threshold} \
            --gene_nt_extend !{params.gene_nt_extend} \
            --unique_threshold !{params.unique_threshold} \
            --threads !{task.cpus} \
            --force \
            --verbose !{noclean} !{spades_options} > ${LOG_DIR}/ariba.out 2> ${LOG_DIR}/ariba.err

    ariba summary !{dataset_name}/summary !{dataset_name}/report.tsv \
            --cluster_cols assembled,match,known_var,pct_id,ctg_cov,novel_var \
            --col_filter n --row_filter n > ${LOG_DIR}/ariba-summary.out 2> ${LOG_DIR}/ariba-summary.err

    rm -rf ariba.tmp*

    if [ "!{params.keep_all_files}" == "false" ]; then
        # Remove Ariba DB that was untarred
        rm -rf !{dataset_name}db
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
    dataset_tarball = path(dataset).getName()
    dataset_name = dataset_tarball.replace('.tar.gz', '')
    """
    mkdir ${dataset_name}
    mkdir ${PROCESS_NAME}
    touch ${dataset_name}/${sample}
    touch ${PROCESS_NAME}/${sample}
    """
}
