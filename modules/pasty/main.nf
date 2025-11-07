process PASTY {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.tsv")        , emit: tsv
    tuple val(meta), path("${prefix}.blastn.tsv") , emit: blast
    tuple val(meta), path("${prefix}.details.tsv"), emit: details
    tuple val(meta), path("*.{log,err}")          , emit: logs, optional: true
    tuple val(meta), path(".command.begin")       , emit: nf_begin
    tuple val(meta), path(".command.err")         , emit: nf_err
    tuple val(meta), path(".command.log")         , emit: nf_log
    tuple val(meta), path(".command.out")         , emit: nf_out
    tuple val(meta), path(".command.run")         , emit: nf_run
    tuple val(meta), path(".command.sh")          , emit: nf_sh
    tuple val(meta), path(".command.trace")       , emit: nf_trace
    tuple val(meta), path("versions.yml")         , emit: versions

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
    """
    pasty \\
        ${task.ext.args} \\
        --prefix $prefix \\
        --input $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pasty: \$(echo \$(pasty --version 2>&1) | sed 's/.*pasty, version //;s/ .*\$//' )
        camlhmp: \$(echo \$(pasty --version 2>&1) | sed 's/.*camlhmp, version //;s/ schema.*\$//' )
    END_VERSIONS
    """
}
