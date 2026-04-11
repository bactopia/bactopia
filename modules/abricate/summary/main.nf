/**
 * Summarize Abricate screening results.
 *
 * Uses [Abricate](https://github.com/tseemann/abricate) to aggregate the
 * per-sample screening reports into a single tab-delimited summary file.
 *
 * @status stable
 * @keywords bacteria, tab-delimited, antimicrobial resistance
 * @tags complexity:simple input-type:single output-type:single features:aggregation
 * @citation abricate
 *
 * @input record(meta, reports)
 * - `meta`: Groovy Record containing aggregation information
 * - `reports`: A collection of Abricate report files from multiple samples
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: Aggregated tab-delimited summary of Abricate results from all samples
 */
nextflow.preview.types = true

process ABRICATE_SUMMARY {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Record,
        reports: Set<Path>
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
        output_dir: task.ext.process_name,
        logs_dir: "${task.ext.process_name}/logs/${task.ext.logs_subdir}/${task.ext.subdir}",
        process_name: task.ext.process_name
    )
    """
    abricate \\
        --summary \\
        ${reports.join(' ')} > ${prefix}.tsv

    # Cleanup

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
    END_VERSIONS
    """
}
