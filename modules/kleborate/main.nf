nextflow.preview.types = true

process KLEBORATE {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fastas) : Tuple<Map, Path>

    output:
    txt      = tuple(meta, file("*.txt"))
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
    mkdir results/
    kleborate \\
        ${task.ext.args} \\
        --outdir results/ \\
        --assemblies ${fastas}

    # Rename output file to include the prefix name
    find results/ -name "*output.txt" -print0 | while read -d \$'\0' file; do mv "\$file" "${prefix}.txt"; done

    # Negative results will not have an output file
    if [ ! -f "${prefix}.txt" ]; then
        touch "${prefix}.txt"
    fi

    # cleanup
    rm -rf results/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kleborate: \$( echo \$(kleborate --version | sed 's/Kleborate v//;'))
    END_VERSIONS
    """
}
