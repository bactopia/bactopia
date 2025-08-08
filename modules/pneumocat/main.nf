process PNEUMOCAT {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}" :
        "${task.ext.docker}" }"

    input:
    tuple val(meta), path(reads)

    when:
    meta.single_end == false

    output:
    tuple val(meta), path("*.xml")                 , optional: true, emit: xml
    tuple val(meta), path("*.coverage_summary.txt"), optional: true, emit: txt
    path "versions.yml"                            , emit: versions
    path ".command.begin"                          , emit: begin
    path ".command.err"                            , emit: err
    path ".command.log"                            , emit: log
    path ".command.out"                            , emit: out
    path ".command.run"                            , emit: run
    path ".command.sh"                             , emit: sh
    path ".command.trace"                          , emit: trace

    script:
    def VERSION = '1.2.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    prefix = task.ext.prefix ?: "${meta.id}"
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
