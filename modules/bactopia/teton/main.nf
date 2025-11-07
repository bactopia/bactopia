process BACTOPIA_SAMPLESHEET {
    tag "${prefix}"
    label 'process_single'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(classification)

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

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/teton/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/teton/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    meta.runtype = _meta.runtype
    meta.teton_reads = _meta.teton_reads
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
        ${task.ext.outdir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sizemeup: \$(echo \$(sizemeup --version 2>&1) | sed 's/.*sizemeup-main, version //;s/ .*\$//' )
    END_VERSIONS
    """
}
