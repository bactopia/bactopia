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
 * @input record(meta, fna)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format
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
    (_meta: Map, fna: Path): Record

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}.tsv"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.tsv")
        ],
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

    def is_compressed = fna.getName().endsWith(".gz") ? true : false
    def fna_name = fna.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fna} > ${fna_name}
    fi

    meningotype \\
        ${task.ext.args} \\
        ${fna_name} \\
        > ${prefix}.tsv

    # Cleanup
    rm -rf ${fna_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meningotype: \$( echo \$(meningotype --version 2>&1) | sed 's/^.*meningotype v//' )
    END_VERSIONS
    """
}
