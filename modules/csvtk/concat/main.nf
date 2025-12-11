/**
 * Concatenate two or more CSV (or TSV) tables into a single table.
 *
 * This process executes csvtk_concat to perform analysis
 *
 * @status stable
 * @keywords concatenate, tsv, csv
 * @tags complexity:simple input-type:multiple output-type:single
 * @citation csvtk_concat
 *
 * @input tuple(meta, csv)
 * - `meta`: Groovy Map containing sample information
 * - `csv`: CSV/TSV formatted files
 *
 * @input in_format
 * Input format (csv, tab, or a delimiting character)
 *
 * @input out_format
 * Output format (csv, tab, or a delimiting character)
 *
 * @output csv      Concatenated CSV/TSV file
 * @output logs     Optional tool execution logs
 * @output nf_logs  Nextflow execution logs
 * @output versions Software version information (YAML format)
 */
nextflow.preview.types = true

process CSVTK_CONCAT {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, csv) : Tuple<Map, Set<Path>>
    in_format    : String
    out_format   : String

    stage:
    stageAs 'inputs/*', csv

    output:
    csv      = tuple(meta, file("${prefix}.${out_extension}"))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, file("versions.yml"))

    script:
    out_extension = out_format == "tsv" ? 'tsv' : 'csv'
    subdir = _meta.subdir ? "${_meta.subdir}/" : ''
    prefix = _meta.id

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
    extra_args = _meta.args ? "${_meta.args}" : ""
    """
    # Create a file of files for csvtk
    ls inputs/ | awk '{ print "inputs/"\$1 }' > fofn.txt

    csvtk \\
        concat \\
        ${task.ext.args} ${extra_args} \\
        --num-cpus ${task.cpus} \\
        ${delimiter}  \\
        ${out_delimiter} \\
        --out-file ${prefix}.${out_extension} \\
        --infile-list fofn.txt

    # Cleanup
    rm -rf fofn.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """
}
