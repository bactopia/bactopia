/**
 * Concatenate two or more CSV (or TSV) tables into a single table.
 *
 * This process executes csvtk_join to perform analysis
 *
 * @status stable
 * @keywords concatenate, tsv, csv
 * @tags complexity:simple input-type:multiple output-type:single
 * @citation csvtk_join
 *
 * @input tuple(meta, csv1, csv2)
 * - `meta`: Groovy Map containing sample information
 * - `csv1`: Input file
 * - `csv2`: Input file
 *
 * @input in_format
 * Input format (csv, tab, or a delimiting character)
 *
 * @input out_format
 * Output format (csv, tab, or a delimiting character)
 *
 * @input key
 * Path parameter for key
 *
 * @output csv      Concatenated CSV/TSV file
 * @output logs     Optional tool execution logs
 * @output nf_logs  Nextflow execution logs
 * @output versions Software version information (YAML format)
 */
nextflow.preview.types = true

process CSVTK_JOIN {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, csv1, csv2) : Tuple<Map, Path, Path>
    in_format           : String
    out_format          : String
    key                 : String

    output:
    csv      = tuple(meta, file("${prefix}.${out_extension}"))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, file("versions.yml"))

    script:
    out_extension = out_format == "tsv" ? 'tsv' : 'csv'
    subdir = _meta.subdir ? "${_meta.subdir}/" : ''
    prefix = task.ext.prefix ?: "${_meta.id}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "merged-results"
    meta.logs_dir = "merged-results/logs/${prefix}-concat/${subdir}"
    meta.process_name = task.ext.process_name
    def delimiter = in_format == "tsv" ? "--tabs" : (in_format == "csv" ? "" : "--delimiter '${in_format}'")
    def out_delimiter = out_format == "tsv" ? "--out-tabs" : (out_format == "csv" ? "" : "--out-delimiter '${out_format}'")
    """
    csvtk \\
        join \\
        ${task.ext.args} \\
        --fields ${key} \\
        --num-cpus ${task.cpus} \\
        ${delimiter}  \\
        ${out_delimiter} \\
        --out-file ${prefix}.${out_extension} \\
        ${csv1} ${csv2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """
}
