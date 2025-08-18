process BACTOPIA_SAMPLESHEET {
    tag "$meta.id"
    label 'process_single'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}${task.ext.singularity_version}" :
        "${task.ext.docker}${task.ext.docker_version}" }"

    input:
    tuple val(meta), path(classification)

    output:
    tuple val(meta), path("${prefix}.bacteria.tsv")   , emit: bacteria_tsv
    tuple val(meta), path("${prefix}.nonbacteria.tsv"), emit: nonbacteria_tsv
    tuple val(meta), path("${prefix}-sizemeup.txt")   , emit: sizemeup
    tuple val(meta), path("*.{log,err}")   , emit: logs, optional: true
    tuple val(meta), path(".command.out")  , emit: nf_out
    tuple val(meta), path(".command.err")  , emit: nf_err
    tuple val(meta), path(".command.log")  , emit: nf_log
    tuple val(meta), path(".command.sh")   , emit: nf_sh
    tuple val(meta), path(".command.trace"), emit: nf_trace
    tuple val(meta), path(".command.run")  , emit: nf_run, optional: true
    tuple val(meta), path(".command.begin"), emit: nf_begin
    tuple val(meta), path("versions.yml")  , emit: versions
    tuple val(meta), path("*-{error,merged}.txt"), optional: true

    script:
    prefix = task.ext.suffix ? "${task.ext.suffix}" : "${meta.id}"
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs"
    meta.process_name = task.ext.process_name
    """
    # determine genome size and create sample sheet
    sizemeup \\
        --query $classification \\
        --prefix ${prefix}

    # create sample sheet
    teton-prepare.py \\
        ${prefix} \\
        ${prefix}-sizemeup.txt \\
        ${meta.runtype} \\
        ${meta.teton_reads} \\
        ${params.outdir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sizemeup: \$(echo \$(sizemeup --version 2>&1) | sed 's/.*sizemeup-main, version //;s/ .*\$//' )
    END_VERSIONS
    """
}
