process SEROBA_RUN {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}" :
        "${task.ext.docker}" }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("results/${prefix}.tsv")                              , emit: tsv
    tuple val(meta), path("results/detailed_serogroup_info.txt"), optional: true, emit: txt
    path "versions.yml"                                                         , emit: versions
    path ".command.begin", emit: begin
    path ".command.err", emit: err
    path ".command.log", emit: log
    path ".command.out", emit: out
    path ".command.run", emit: run
    path ".command.sh", emit: sh
    path ".command.trace", emit: trace

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    seroba \\
        runSerotyping \\
        $reads \\
        results $args

    # Avoid name collisions
    mv results/pred.tsv results/${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seroba: \$(seroba version)
    END_VERSIONS
    """
}
