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
    logs     = tuple(meta, file("*.{log,err}", optional: true))
    nf_begin = tuple(meta, file(".command.begin"))
    nf_err   = tuple(meta, file(".command.err"))
    nf_log   = tuple(meta, file(".command.log"))
    nf_out   = tuple(meta, file(".command.out"))
    nf_run   = tuple(meta, file(".command.run"))
    nf_sh    = tuple(meta, file(".command.sh"))
    nf_trace = tuple(meta, file(".command.trace"))
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
