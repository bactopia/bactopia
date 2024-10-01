// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'ariba')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::ariba=2.14.6"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process ARIBA_RUN {
    tag "$meta.id - $db_name"
    label 'process_low'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ariba:2.14.6--py39h67e14b5_4' :
        'quay.io/biocontainers/ariba:2.14.6--py39h67e14b5_4' }"

    input:
    tuple val(meta), path(reads)
    path(db)

    output:
    tuple val(meta), path("results/*")                    , emit: results
    tuple val(meta), path("results/${prefix}-report.tsv") , emit: report
    tuple val(meta), path("results/${prefix}-summary.csv"), emit: summary
    path "*.{log,err}"                                    , emit: logs, optional: true
    path ".command.*"                                     , emit: nf_logs
    path "versions.yml"                                   , emit: versions

    when:
    meta.single_end == false

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    db_name = db.getName().replace('.tar.gz', '')
    """
    tar -xzvf ${db}
    mv ${db_name} ${db_name}db
    ariba \\
        run \\
        ${db_name}db/ \\
        ${reads} \\
        ${db_name} \\
        $options.args \\
        --threads $task.cpus

    ariba \\
        summary \\
        ${db_name}/summary \\
        ${db_name}/report.tsv \\
        --cluster_cols assembled,match,known_var,pct_id,ctg_cov,novel_var \\
        --col_filter n \\
        --row_filter n

    # Rename to avoid naming collisions
    mv ${db_name}/report.tsv ${db_name}/${prefix}-report.tsv
    mv ${db_name}/summary.csv ${db_name}/${prefix}-summary.csv
    mv ${db_name}/ results/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ariba:  \$(echo \$(ariba version 2>&1) | sed 's/^.*ARIBA version: //;s/ .*\$//')
    END_VERSIONS
    """
}
