/**
 * Serotyping and finetyping of *Neisseria meningitidis*.
 *
 * Uses [Meningotype](https://github.com/MDU-PHL/meningotype) to predict the serogroup (capsule),
 * PorA variable regions, and FetA variable regions of *N. meningitidis* assemblies. This provides
 * a comprehensive molecular typing profile (e.g., B:P1.7-2,4:F1-5) used in surveillance.
 *
 * @status stable
 * @keywords bacteria, neisseria meningitidis, serotype, serogroup, finetyping, pora, feta, capsule
 * @tags complexity:simple input-type:single output-type:single features:conditional-logic
 * @citation meningotype
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: Tab-delimited meningotype results with serogroup, PorA, and FetA predictions
 */
nextflow.preview.types = true

process MENINGOTYPE {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, assembly: Path): Record

    output:
    record(
        meta: meta,
        tsv: file("${prefix}.tsv"),
        results: [file("${prefix}.tsv")],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

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

    def is_compressed = assembly.getName().endsWith(".gz") ? true : false
    def assembly_name = assembly.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${assembly} > ${assembly_name}
    fi

    meningotype \\
        ${task.ext.args} \\
        ${assembly_name} \\
        > ${prefix}.tsv

    # Cleanup
    rm -rf ${assembly_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meningotype: \$( echo \$(meningotype --version 2>&1) | sed 's/^.*meningotype v//' )
    END_VERSIONS
    """
}
