process PNEUMOCAT {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.xml")                      , emit: xml, optional: true
    tuple val(meta), path("*.coverage_summary.txt")     , emit: txt, optional: true
    tuple val(meta), path("*.{log,err}")                , emit: logs, optional: true
    tuple val(meta), path(".command.begin")             , emit: nf_begin
    tuple val(meta), path(".command.err")               , emit: nf_err
    tuple val(meta), path(".command.log")               , emit: nf_log
    tuple val(meta), path(".command.out")               , emit: nf_out
    tuple val(meta), path(".command.run")               , emit: nf_run
    tuple val(meta), path(".command.sh")                , emit: nf_sh
    tuple val(meta), path(".command.trace")             , emit: nf_trace
    tuple val(meta), path("versions.yml")               , emit: versions

    script:
    def VERSION = '1.2.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    prefix = task.ext.prefix ?: "${meta.id}"
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs"
    meta.process_name = task.ext.process_name
    """
    PneumoCaT.py \\
        --input_directory ./ \\
        --threads $task.cpus \\
        --output_dir ./

    # clean up
    rm -rf *.bam *.bai ComponentComplete.txt

    # PneumoCAT uses first match in a glob, so moves between R1 and R2
    if [ -f ${prefix}_R1.results.xml ]; then
        mv ${prefix}_R1.results.xml ${prefix}.results.xml
    else
        mv ${prefix}_R2.results.xml ${prefix}.results.xml
    fi
    mv logs/* ./
    mv coverage_summary.txt ${prefix}.coverage_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pneumocat: $VERSION
    END_VERSIONS
    """
}
