nextflow.preview.types = true

process BTYPER3 {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta) : Tuple<Map, Path>

    output:
    tsv          = tuple(meta, file("${prefix}.tsv"))
    supplemental = tuple(meta, file("supplemental/*", optional: true))
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
    # Btyper3 does not accept compressed files
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    btyper3 \\
        ${task.ext.args} \\
        --output ./ \\
        --input ${fasta_name}

    mv btyper3_final_results/ supplemental/
    mv supplemental/${prefix}_final_results.txt ./${prefix}.tsv

    # Cleanup
    rm -rf ${fasta_name} ${fasta_name}.njs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        btyper3: \$(echo \$(btyper3 --version 2>&1) | sed 's/^.*btyper3 //;' ))
    END_VERSIONS
    """
}
