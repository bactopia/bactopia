process CSVTK_CONCAT {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(csv, stageAs: 'inputs/*')
    val in_format
    val out_format

    output:
    tuple val(meta), path("${prefix}.${out_extension}"), emit: csv
    tuple val(meta), path("*.{log,err}")               , emit: logs, optional: true
    tuple val(meta), path(".command.begin")            , emit: nf_begin
    tuple val(meta), path(".command.err")              , emit: nf_err
    tuple val(meta), path(".command.log")              , emit: nf_log
    tuple val(meta), path(".command.out")              , emit: nf_out
    tuple val(meta), path(".command.run")              , emit: nf_run
    tuple val(meta), path(".command.sh")               , emit: nf_sh
    tuple val(meta), path(".command.trace")            , emit: nf_trace
    tuple val(meta), path("versions.yml")              , emit: versions

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
        $task.ext.args $extra_args \\
        --num-cpus $task.cpus \\
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
