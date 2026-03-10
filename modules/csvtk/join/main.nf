/**
 * Join two CSV or TSV files based on common fields.
 *
 * Uses [csvtk join](https://github.com/shenwei356/csvtk) to merge two tabular files horizontally
 * by matching values in a specified key column (similar to a SQL JOIN). It supports inner, left,
 * right, and outer joins via optional arguments.
 *
 * @status stable
 * @keywords utility, table, join, merge, csv, tsv, csvtk, relational
 * @tags complexity:simple input-type:multiple output-type:single features:conditional-logic
 * @citation csvtk
 *
 * @input record(_meta, csv1, csv2)
 * - `_meta`: Groovy Map containing sample information
 * - `csv1`: The first CSV/TSV file (Left table)
 * - `csv2`: The second CSV/TSV file (Right table)
 *
 * @input in_format
 * Input format string ('csv', 'tsv', or a specific delimiter character)
 *
 * @input out_format
 * Output format string ('csv', 'tsv', or a specific delimiter character)
 *
 * @input key
 * The column name(s) or index(es) to use as the join key (e.g., "sample_id" or "1")
 *
 * @output record(meta, csv, results, logs, nf_logs, versions)
 * - `csv`: The joined tabular file (*.csv or *.tsv)
 */
nextflow.preview.types = true

process CSVTK_JOIN {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta: Map, csv1: Path, csv2: Path): Record
    in_format : String
    out_format: String
    key       : String

    output:
    record(
        meta: meta,
        csv: file("${prefix}.${out_extension}"),
        results: [file("${prefix}.${out_extension}")],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

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
