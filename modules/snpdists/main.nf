nextflow.preview.types = true

process SNPDISTS {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, alignment) : Tuple<Map, Path>

    output:
    tsv      = tuple(meta, file("*.tsv"))
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
    process_name = _meta.process_name ?: task.ext.process_name

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.ext.process_name}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = ""
    meta.logs_dir = "${process_name}/logs/"
    meta.process_name = process_name
    """
    snp-dists \\
        ${task.ext.args} \\
        ${alignment} > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snp-dists: \$(snp-dists -v 2>&1 | sed 's/snp-dists //;')
    END_VERSIONS
    """
}
