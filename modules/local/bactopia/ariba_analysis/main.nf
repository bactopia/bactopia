nextflow.enable.dsl = 2

// Assess cpu and memory of current system
include { get_resources; save_files } from '../../utilities/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)
PROCESS_NAME = "ariba_analysis"

process ARIBA_ANALYSIS {
    /* Run reads against all available (if any) ARIBA datasets */
    tag "${sample} - ${dataset_name}"
    label PROCESS_NAME

    publishDir "${params.outdir}/${sample}",
        mode: params.publish_mode,
        overwrite: params.force,
        saveAs: { filename -> save_files(filename:filename, process_name:PROCESS_NAME, logs_subdir:dataset_name) }

    input:
    tuple val(sample), val(single_end), path(fq)
    each path(dataset)

    output:
    path "${dataset_name}/*", emit: results
    path "*.std{out,err}.txt", emit: logs
    path ".command.*", emit: nf_logs
    path "*.version.txt", emit: version

    when:
    single_end == false

    shell:
    dataset_tarball = dataset.getName()
    dataset_name = dataset_tarball.replace('.tar.gz', '')
    spades_options = params.spades_options ? "--spades_options '${params.spades_options}'" : ""
    noclean = params.ariba_no_clean ? "--noclean" : ""
    '''
    tar -xzvf !{dataset_tarball}
    mv !{dataset_name} !{dataset_name}db
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
            --verbose !{noclean} !{spades_options} > ariba.stdout.txt 2> ariba.stderr.txt

    ariba summary !{dataset_name}/summary !{dataset_name}/report.tsv \
            --cluster_cols assembled,match,known_var,pct_id,ctg_cov,novel_var \
            --col_filter n --row_filter n > ariba-summary.stdout.txt 2> ariba-summary.stderr.txt

    # Cleanup
    rm -rf ariba.tmp*
    rm -rf !{dataset_name}db

    # Capture version
    ariba version >> ariba.version.txt 2>&1
    '''

    stub:
    dataset_tarball = path(dataset).getName()
    dataset_name = dataset_tarball.replace('.tar.gz', '')
    """
    mkdir ${dataset_name}
    touch ${dataset_name}/${sample}
    """
}
