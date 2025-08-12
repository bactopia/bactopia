process CSVTK_CONCAT {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(csv, stageAs: 'inputs/*')
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
    prefix = task.ext.prefix ? "${task.ext.prefix}" : "${meta.id}"
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs"
    meta.process_name = task.ext.process_name
    def delimiter = in_format == "tsv" ? "--tabs" : (in_format == "csv" ? "" : "--delimiter '${in_format}'")
    def out_delimiter = out_format == "tsv" ? "--out-tabs" : (out_format == "csv" ? "" : "--out-delimiter '${out_format}'")
    out_extension = out_format == "tsv" ? 'tsv' : 'csv'
    """
    # Create a file of files for csvtk
    ls inputs/ | awk '{ print "inputs/"\$1 }' > fofn.txt

    csvtk \\
        concat \\
        $task.ext.args \\
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
