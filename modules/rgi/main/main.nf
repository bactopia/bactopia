process RGI_MAIN {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}" :
        "${task.ext.docker}" }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.json"), emit: json, optional: true
    tuple val(meta), path("*.txt") , emit: tsv
    path "versions.yml"            , emit: versions
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
    rgi \\
        main \\
        $args \\
        --clean \\
        --data wgs \\
        --num_threads $task.cpus \\
        --output_file $prefix \\
        --input_sequence $fasta

    # Remove empty json files
    if grep "^{}\$" ${prefix}.json; then
        rm ${prefix}.json
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(rgi main --version)
    END_VERSIONS
    """
}
