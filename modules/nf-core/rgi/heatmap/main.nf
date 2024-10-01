// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'rgi_heatmap')
options.btype = options.btype ?: "comparative"
conda_tools   = "bioconda::rgi=6.0.3"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process RGI_HEATMAP {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rgi:6.0.3--pyha8f3691_0' :
        'quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0' }"

    input:
    tuple val(meta), path(json, stageAs: 'json/*')

    output:
    tuple val(meta), path("*.{csv,eps,png}"), emit: heatmap, optional: true
    path "*.{log,err}"                      , emit: logs, optional: true
    path ".command.*"                       , emit: nf_logs
    path "versions.yml"                     , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    NUM_SAMPLES=\$(ls json/ | wc -l)
    if [[ "\${NUM_SAMPLES}" -gt 1 ]]; then
        rgi \\
            heatmap \\
            $options.args \\
            --output $prefix \\
            --input json/
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(rgi main --version)
    END_VERSIONS
    """
}
