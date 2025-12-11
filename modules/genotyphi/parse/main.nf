/**
 * Parse GenoTyphi results from Mykrobe JSON output.
 *
 * This process executes genotyphi_parse to perform analysis
 *
 * @status stable
 * @keywords salmonella, typhi, genotyping, parse
 * @tags complexity:simple input-type:single output-type:single
 * @citation genotyphi_parse
 *
 * @input tuple(meta, json)
 * - `meta`: Groovy Map containing sample information
 * - `json`: Mykrobe JSON output file
 *
 * @output tsv      Tab-delimited genotyphi results
 * @output logs     Optional tool execution logs
 * @output nf_logs  Nextflow execution logs
 * @output versions Software version information (YAML format)
 */
nextflow.preview.types = true

process GENOTYPHI_PARSE {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, json) : Tuple<Map, Path>

    output:
    tsv      = tuple(meta, files("*.tsv"))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, file("versions.yml"))

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
    """
    parse_typhi_mykrobe.py \\
        --jsons ${json} \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genotyphi: \$(echo \$(genotyphi --version 2>&1) | sed 's/^.*GenoTyphi v//;' )
    END_VERSIONS
    """
}
