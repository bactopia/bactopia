process LEGSTA {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}":"${task.ext.docker}" }"

    input:
    tuple val(meta), path(seqs)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "*.{log,err}"            , emit: logs, optional: true
    path ".command.{begin,err,log,out,run,sh,trace}", emit: nf_logs
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "task.ext.args: ${task.ext.args}"
    
    # Create .command.begin
    date > .command.begin
    
    legsta \\
        $args \\
        $seqs | sed 's/.fna//; s/.gz//' > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        legsta: \$(echo \$(legsta --version 2>&1) | sed 's/^.*legsta //; s/ .*\$//;')
    END_VERSIONS
    """
}
