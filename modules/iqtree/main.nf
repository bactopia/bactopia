nextflow.preview.types = true

process IQTREE {
    tag "${prefix}"
    label 'process_medium'
    label 'process_long'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, alignment) : Tuple<Map, Path>

    output:
    supplemental = tuple(meta, file("${process_name}/*"))
    phylogeny    = tuple(meta, file(treefile))
    alignment    = tuple(meta, alignment)
    aln_tree     = tuple(meta, alignment, file(treefile))
    logs         = tuple(meta, file("*.{log,err}", optional: true))
    nf_begin     = tuple(meta, file(".command.begin"))
    nf_err       = tuple(meta, file(".command.err"))
    nf_log       = tuple(meta, file(".command.log"))
    nf_out       = tuple(meta, file(".command.out"))
    nf_run       = tuple(meta, file(".command.run"))
    nf_sh        = tuple(meta, file(".command.sh"))
    nf_trace     = tuple(meta, file(".command.trace"))
    versions     = tuple(meta, file("versions.yml"))

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"
    process_name = _meta.process_name == "iqtree-fast" ? "iqtree-fast" : task.ext.process_name
    args = process_name == "iqtree-fast" ? task.ext.fast_args : task.ext.args
    treefile = process_name == "iqtree-fast" ? "${process_name}/${prefix}.treefile" : "${prefix}.treefile"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = ""
    meta.logs_dir = "${process_name}/logs/"
    meta.process_name = process_name
    """
    iqtree \\
        ${args} \\
        -s ${alignment} \\
        -nt ${task.cpus} \\
        -ntmax ${task.cpus} \\
        -pre ${prefix}

    # Only gzip files if they exist
    if [[ -f "${prefix}.alninfo" ]]; then
        gzip ${prefix}.alninfo
    fi

    mkdir temp
    mv ${prefix}* temp/
    mv temp/ ${process_name}/

    if [ "${process_name}" != "iqtree-fast" ]; then
        mv ${process_name}/${prefix}.treefile ./
        mv ${process_name}/${alignment} ./
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
    END_VERSIONS
    """
}
