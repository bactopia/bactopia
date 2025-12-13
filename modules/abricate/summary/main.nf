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
 * @input tuple(meta, reports)
 * - `meta`: Groovy Map containing sample information
 * - `reports`: A collection of Abricate report files from multiple samples
 *
 * @output report   A merged tab-delimited file with [Abricate](https://github.com/tseemann/abricate) results from all samples
 * @output logs     Optional software execution logs containing warnings/errors
 * @output nf_logs  Nextflow execution scripts and logs for debugging
 * @output versions A YAML formatted file with software versions
 */
nextflow.preview.types = true

process ABRICATE_SUMMARY {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, reports): Tuple<Map, Set<Path>>

    output:
    report   = tuple(meta, file("${prefix}.tsv"))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, files("versions.yml"))

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
