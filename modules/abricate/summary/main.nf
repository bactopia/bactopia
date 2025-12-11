/**
 * Screen assemblies for antimicrobial resistance against multiple databases.
 *
 * This process executes abricate_summary to perform analysis
 *
 * @status stable
 * @keywords bacteria, assembly, antimicrobial reistance
 * @tags complexity:simple input-type:single output-type:single
 * @citation abricate_summary
 *
 * @input tuple(meta, reports)
 * - `meta`: Groovy Map containing sample information
 * - `reports`: Input file
 *
 * @output report   Report
 * @output logs     Optional tool execution logs
 * @output nf_logs  Nextflow execution logs
 * @output versions Software version information (YAML format)
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
    versions = tuple(meta, file("versions.yml"))

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
