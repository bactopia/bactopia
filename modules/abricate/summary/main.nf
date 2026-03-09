/**
 * Summarize Abricate screening results.
 *
 * Uses [Abricate](https://github.com/tseemann/abricate) to aggregate the
 * per-sample screening reports into a single tab-delimited summary file.
 *
 * @status stable
 * @keywords bacteria, tab-delimited, antimicrobial resistance
 * @tags complexity:simple input-type:single output-type:single
 * @citation abricate
 *
 * @input record(meta, reports)
 * - `meta`: Groovy Map containing aggregation information
 * - `reports`: A collection of Abricate report files from multiple samples
 *
 * @output record(meta, report, results, logs, nf_logs, versions)
 */
nextflow.preview.types = true

process ABRICATE_SUMMARY {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, reports: Set<Path>): Record

    output:
    record(
        meta: meta,
        report: file("${prefix}.tsv"),
        results: [file("${prefix}.tsv")],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    prefix = task.ext.prefix ?: "${_meta.id}"
    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${task.ext.process_name}"
    meta.logs_dir = "${task.ext.process_name}/logs/${task.ext.logs_subdir}/${task.ext.subdir}"
    meta.process_name = task.ext.process_name
    """
    abricate \\
        --summary \\
        ${reports} > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
    END_VERSIONS
    """
}
