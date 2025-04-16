process ABRICATE_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda task.ext.env.conda
    container task.ext.env.container

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.txt"), emit: report
    path "*.{log,err}"            , emit: logs, optional: true
    path ".command.*"             , emit: nf_logs
    path "versions.yml"           , emit: versions

    script:
    prefix = task.ext.prefix ? "${task.ext.prefix}" : "${meta.id}"
    """
    abricate \\
        $assembly \\
        $task.ext.args \\
        --threads $task.cpus > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
    END_VERSIONS
    """
}
