process SHIGATYPER {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}" :
        "${task.ext.docker}" }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.tsv")     , emit: tsv
    tuple val(meta), path("${prefix}-hits.tsv"), emit: hits, optional: true
    path "versions.yml"                        , emit: versions
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

    if (meta.runtype == "ont") {
        """
        shigatyper \\
            $args \\
            --SE $reads \\
            --ont \\
            --name $prefix

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            shigatyper: \$(echo \$(shigatyper --version 2>&1) | sed 's/^.*ShigaTyper //' )
        END_VERSIONS
        """
    } else if (meta.single_end) {
        """
        shigatyper \\
            $args \\
            --SE $reads \\
            --name $prefix

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            shigatyper: \$(echo \$(shigatyper --version 2>&1) | sed 's/^.*ShigaTyper //' )
        END_VERSIONS
        """
    } else {
        """
        shigatyper \\
            $args \\
            --R1 ${reads[0]} \\
            --R2 ${reads[1]} \\
            --name $prefix

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            shigatyper: \$(echo \$(shigatyper --version 2>&1) | sed 's/^.*ShigaTyper //' )
        END_VERSIONS
        """
    }
}
