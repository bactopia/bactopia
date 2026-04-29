/**
 * Run all available MLST schemes for a species against an assembly
 *
 * Uses [GigaTyper](https://github.com/rpetit3/gigatyper) to run all available mlst schemes for a species against an assembly.
 *
 * @status stable
 * @keywords mlst, typing, multi-scheme
 * @tags complexity:simple input-type:single output-type:single features:
 * @citation gigatyper
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Record containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: MLST results across all schemes
 */
nextflow.enable.types = true

process GIGATYPER {
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

    // Create a new meta record
    meta = record(
        id: "${prefix}-${task.process}",
        name: prefix,
        scope: task.ext.scope,
        output_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}",
        logs_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}",
        process_name: task.ext.process_name
    )

    """
    gigatyper \\
        --input ${fna} \\
        --prefix ${prefix} \\
        --threads ${task.cpus} \\
        ${task.ext.args} \\
        > ${prefix}.tsv

    # Cleanup

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gigatyper: \$( gigatyper --version 2>&1 | sed 's/gigatyper //' )
    END_VERSIONS
    """
}
