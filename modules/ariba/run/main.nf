process ARIBA_RUN {
    tag "$meta.id - $db_name"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(reads)
    path(db)

    output:
    tuple val(meta), path("results/*")                    , emit: results
    tuple val(meta), path("results/${prefix}-report.tsv") , emit: report
    tuple val(meta), path("results/${prefix}-summary.csv"), emit: summary
    tuple val(meta), path("*.{log,err}")                  , emit: logs, optional: true
    tuple val(meta), path(".command.begin")               , emit: nf_begin
    tuple val(meta), path(".command.err")                 , emit: nf_err
    tuple val(meta), path(".command.log")                 , emit: nf_log
    tuple val(meta), path(".command.out")                 , emit: nf_out
    tuple val(meta), path(".command.run")                 , emit: nf_run
    tuple val(meta), path(".command.sh")                  , emit: nf_sh
    tuple val(meta), path(".command.trace")               , emit: nf_trace
    tuple val(meta), path("versions.yml")                 , emit: versions

    script:
    prefix = task.ext.prefix ? "${task.ext.prefix}" : "${meta.id}"
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs"
    meta.process_name = task.ext.process_name
    db_name = db.getName().replace('.tar.gz', '')
    """
    tar -xzvf ${db}
    mv ${db_name} ${db_name}db
    ariba \\
        run \\
        ${db_name}db/ \\
        ${reads} \\
        ${db_name} \\
        $task.ext.args \\
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

    # Cleanup
    rm -rf ${db_name}db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ariba:  \$(echo \$(ariba version 2>&1) | sed 's/^.*ARIBA version: //;s/ .*\$//')
    END_VERSIONS
    """
}
