process RGI_HEATMAP {
    tag "$meta.id"
    label 'process_single'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}" :
        "${task.ext.docker}" }"

    input:
    tuple val(meta), path(json, stageAs: 'json/*')

    output:
    tuple val(meta), path("*.{csv,eps,png}"), emit: heatmap, optional: true
    path "versions.yml"                     , emit: versions
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
    NUM_SAMPLES=\$(ls json/ | wc -l)
    if [[ "\${NUM_SAMPLES}" -gt 1 ]]; then
        rgi \\
            heatmap \\
            $args \\
            --output $prefix \\
            --input json/
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(rgi main --version)
    END_VERSIONS
    """
}
