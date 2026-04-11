/**
 * In silico Sequence Based Typing (SBT) of *Legionella pneumophila*.
 *
 * Uses [Legsta](https://github.com/tseemann/legsta) to determine the Sequence Based Type (SBT)
 * of *L. pneumophila* isolates. It aligns the assembly against the standard 7-gene schema
 * (flaA, pilE, asd, mip, mompS, proA, neuA) to assign allele numbers and the resulting Sequence Type.
 *
 * @status stable
 * @keywords bacteria, legionella, pneumophila, typing, sbt, mlst, serogroup
 * @tags complexity:simple input-type:single output-type:single features:database-dependent
 * @citation legsta
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Record containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: Tab-delimited Legionella pneumophila SBT results with allele numbers and sequence type
 */
nextflow.preview.types = true

process LEGSTA {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Record,
        fna: Path
    )

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
    def _meta = meta
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = record(
        id: "${prefix}-${task.process}",
        name: prefix,
        scope: task.ext.scope,
        output_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}",
        logs_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}",
        process_name: task.ext.process_name
    )
    """
    legsta \\
        ${task.ext.args} \\
        ${fna} | sed 's/.fna//; s/.gz//' > ${prefix}.tsv

    # Cleanup

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        legsta: \$(echo \$(legsta --version 2>&1) | sed 's/^.*legsta //; s/ .*\$//;')
    END_VERSIONS
    """
}
