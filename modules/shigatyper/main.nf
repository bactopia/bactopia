process SHIGATYPER {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.tsv")     , emit: tsv
    tuple val(meta), path("${prefix}-hits.tsv"), emit: hits, optional: true
    tuple val(meta), path("*.{log,err}")       , emit: logs, optional: true
    tuple val(meta), path(".command.begin")    , emit: nf_begin
    tuple val(meta), path(".command.err")      , emit: nf_err
    tuple val(meta), path(".command.log")      , emit: nf_log
    tuple val(meta), path(".command.out")      , emit: nf_out
    tuple val(meta), path(".command.run")      , emit: nf_run
    tuple val(meta), path(".command.sh")       , emit: nf_sh
    tuple val(meta), path(".command.trace")    , emit: nf_trace
    tuple val(meta), path("versions.yml")      , emit: versions

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name

    if (_meta.runtype == "ont") {
        """
        shigatyper \\
            ${task.ext.args} \\
            --SE $reads \\
            --ont \\
            --name $prefix

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            shigatyper: \$(echo \$(shigatyper --version 2>&1) | sed 's/^.*ShigaTyper //' )
        END_VERSIONS
        """
    } else if (_meta.single_end) {
        """
        shigatyper \\
            ${task.ext.args}  \\
            --SE $reads \\
            --name $prefix

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            shigatyper: \$(echo \$(shigatyper --version 2>&1) | sed 's/^.*ShigaTyper //' )
        END_VERSIONS
        """
    } else {
        """
        shigatyper \\
            ${task.ext.args}  \\
            --R1 ${reads[0]} \\
            --R2 ${reads[1]} \\
            --name $prefix

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            shigatyper: \$(echo \$(shigatyper --version 2>&1) | sed 's/^.*ShigaTyper //' )
        END_VERSIONS
        """
    }
}
