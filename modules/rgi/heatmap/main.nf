process RGI_HEATMAP {
    tag "${prefix}"
    label 'process_single'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(json, stageAs: 'json/*')

    output:
    tuple val(meta), path("*.{csv,eps,png}"), emit: heatmap, optional: true
    tuple val(meta), path("*.{log,err}")    , emit: logs, optional: true
    tuple val(meta), path(".command.begin") , emit: nf_begin
    tuple val(meta), path(".command.err")   , emit: nf_err
    tuple val(meta), path(".command.log")   , emit: nf_log
    tuple val(meta), path(".command.out")   , emit: nf_out
    tuple val(meta), path(".command.run")   , emit: nf_run
    tuple val(meta), path(".command.sh")    , emit: nf_sh
    tuple val(meta), path(".command.trace") , emit: nf_trace
    tuple val(meta), path("versions.yml")   , emit: versions

    script:
    prefix = task.ext.prefix ?: "${_meta.id}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "merged-results/"
    meta.logs_dir = "merged-results/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    """
    NUM_SAMPLES=\$(ls json/ | wc -l)
    if [[ "\${NUM_SAMPLES}" -gt 1 ]]; then
        rgi \\
            heatmap \\
            ${task.ext.args} \\
            --output $prefix \\
            --input json/
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(rgi main --version)
    END_VERSIONS
    """
}
