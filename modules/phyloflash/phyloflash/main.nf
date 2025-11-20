nextflow.preview.types = true

process PHYLOFLASH {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, reads) : Tuple<Map, Path>
    _silva_db       : Path
    _univec_db      : Path

    output:
    supplemental = tuple(meta, file("${prefix}/*"))
    aln          = file("${prefix}/${prefix}.toalign.fasta", optional: true)
    summary      = file("${prefix}/${prefix}.phyloFlash.json", optional: true)
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

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    def read_opts = meta.single_end ? "-read1 ${reads[0]}" : "-read1 ${reads[0]} -read2 ${reads[1]}"
    """
    mkdir ${prefix}
    phyloFlash.pl \\
        ${task.ext.args} \\
        ${read_opts} \\
        -lib ${prefix} \\
        -dbhome . \\
        -CPUs ${task.cpus}

    jsonify-phyloflash.py ${prefix}.phyloFlash > ${prefix}.phyloFlash.json
    mv ${prefix}.* ${prefix}


    if phyloflash-summary.py ${prefix}/ | grep -q -c "WARNING: Multiple SSUs were assembled by SPAdes"; then
        MULTI="1"
    fi

    if [ "${task.ext.allow_multiple_16s}" == "true" ]; then
        MULTI="0"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phyloFlash: \$(echo \$(phyloFlash.pl -version 2>&1) | sed "s/^.*phyloFlash v//")
    END_VERSIONS
    """
}
