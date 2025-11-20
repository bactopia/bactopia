nextflow.preview.types = true

process ARIBA_RUN {
    tag "${prefix} - ${db_name}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, reads) : Tuple<Map, Path>
    db             : Path

    output:
    report       = tuple(meta, file("${prefix}-report.tsv"))
    summary      = tuple(meta, file("${prefix}-summary.csv"))
    supplemental = tuple(meta, file("supplemental/*"))
    logs         = tuple(meta, file("*.{log,err}", optional: true))
    nf_begin     = tuple(meta, file(".command.begin"))
    nf_err       = tuple(meta, file(".command.err"))
    nf_log       = tuple(meta, file(".command.log"))
    nf_out       = tuple(meta, file(".command.out"))
    nf_run       = tuple(meta, file(".command.run"))
    nf_sh        = tuple(meta, file(".command.sh"))
    nf_trace     = tuple(meta, file(".command.trace"))
    versions     = tuple(meta, file("versions.yml"))

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
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
        ${task.ext.args} \\
        --threads ${task.cpus}

    ariba \\
        summary \\
        ${db_name}/summary \\
        ${db_name}/report.tsv \\
        --cluster_cols assembled,match,known_var,pct_id,ctg_cov,novel_var \\
        --col_filter n \\
        --row_filter n

    # Rename to avoid naming collisions
    mv ${db_name}/report.tsv ./${prefix}-report.tsv
    mv ${db_name}/summary.csv ./${prefix}-summary.csv
    mv ${db_name}/ supplemental/

    # Cleanup
    rm -rf ${db_name}db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ariba:  \$(echo \$(ariba version 2>&1) | sed 's/^.*ARIBA version: //;s/ .*\$//')
    END_VERSIONS
    """
}
