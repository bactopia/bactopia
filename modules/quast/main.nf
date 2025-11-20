nextflow.preview.types = true

process QUAST {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta, meta_file) : Tuple<Map, Path, Path>

    output:
    tsv          = tuple(meta, file("${prefix}.tsv"))
    supplemental = tuple(meta, file("supplemental/*"))
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
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    est_ref_size=""
    # Use rev to get the last column easily, then re-reverse it
    ref_size=\$(tail -n 1 ${meta_file} | rev | cut -f 1 | rev)
    if [ "\${ref_size}" != "0" ]; then
        est_ref_size="--est-ref-size \${ref_size}"
    fi

    quast ${fasta_name} \${est_ref_size} \\
        -o supplemental \\
        --threads ${task.cpus} \\
        ${task.ext.args} \\
        --glimmer

    mv supplemental/quast.log ./
    mv supplemental/transposed_report.tsv ${prefix}.tsv

    # Cleanup
    rm -rf ${fasta_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
