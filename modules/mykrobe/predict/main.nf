nextflow.preview.types = true

process MYKROBE_PREDICT {
    tag "${meta.name}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, seqs) : Tuple<Map, Path>
    species       : String

    output:
    csv      = tuple(meta, file("${prefix}.csv"))
    json     = tuple(meta, file("${prefix}.json"))
    logs     = tuple(meta, file("*.{log,err}", optional: true))
    nf_begin = tuple(meta, file(".command.begin"))
    nf_err   = tuple(meta, file(".command.err"))
    nf_log   = tuple(meta, file(".command.log"))
    nf_out   = tuple(meta, file(".command.out"))
    nf_run   = tuple(meta, file(".command.run"))
    nf_sh    = tuple(meta, file(".command.sh"))
    nf_trace = tuple(meta, file(".command.trace"))
    versions = tuple(meta, file("versions.yml"))

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
    is_ont = meta.runtype == "ont" ? "--ont" : ""
    """
    mykrobe \\
        predict \\
        ${task.ext.args} ${is_ont} \\
        --species ${species} \\
        --threads ${task.cpus} \\
        --sample ${prefix} \\
        --format json_and_csv \\
        --output ${prefix} \\
        --seq ${seqs}

    # Cleanup
    rm -rf mykrobe

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mykrobe: \$(echo \$(mykrobe --version 2>&1) | sed 's/^.*mykrobe v//' )
    END_VERSIONS
    """
}
