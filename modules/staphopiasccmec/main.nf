/**
 * Primer based SCCmec typing of S. aureus genomes.
 *
 * Uses [Staphopia SCCmec](https://github.com/staphopia/staphopia-sccmec) to determine the
 * SCCmec type of *Staphylococcus aureus* assemblies. It includes a primer-based approach
 * to identify SCCmec types I-XI.
 *
 * @status stable
 * @keywords staphylococcus aureus, sccmec, typing, mrsa, primers
 * @tags complexity:simple input-type:single output-type:single features:compression,conditional-logic
 * @citation staphopiasccmec
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: TSV file with SCCmec typing results
 */
nextflow.preview.types = true

process STAPHOPIASCCMEC {
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

    staphopia-sccmec --assembly ${fna_name} ${task.ext.args} > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        staphopia-sccmec: \$(staphopia-sccmec --version 2>&1 | sed 's/^.*staphopia-sccmec //')
    END_VERSIONS
    """
}
