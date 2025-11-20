nextflow.preview.types = true

process BACTOPIA_SAMPLESHEET {
    tag "${prefix}"
    label 'process_single'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, classification) : Tuple<Map, Path>

    output:
    bacteria_tsv    = tuple(meta, file("${prefix}.bacteria.tsv"))
    nonbacteria_tsv = tuple(meta, file("${prefix}.nonbacteria.tsv"))
    sizemeup        = tuple(meta, file("${prefix}-sizemeup.txt"))
    logs            = tuple(meta, file("*.{log,err}", optional: true))
    nf_out          = tuple(meta, file(".command.out"))
    nf_err          = tuple(meta, file(".command.err"))
    nf_log          = tuple(meta, file(".command.log"))
    nf_sh           = tuple(meta, file(".command.sh"))
    nf_trace        = tuple(meta, file(".command.trace"))
    nf_run          = tuple(meta, file(".command.run", optional: true))
    nf_begin        = tuple(meta, file(".command.begin"))
    versions        = tuple(meta, file("versions.yml"))

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
        --query ${classification} \\
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
