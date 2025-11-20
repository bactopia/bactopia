nextflow.preview.types = true

process RGI_MAIN {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta) : Tuple<Map, Path>

    output:
    json     = tuple(meta, file("*.json", optional: true))
    tsv      = tuple(meta, file("*.txt"))
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
    """
    rgi \\
        main \\
        ${task.ext.args} \\
        --clean \\
        --data wgs \\
        --num_threads ${task.cpus} \\
        --output_file ${prefix} \\
        --input_sequence ${fasta}

    # Remove empty json files
    if grep "^{}\$" ${prefix}.json; then
        rm ${prefix}.json
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(rgi main --version)
    END_VERSIONS
    """
}
