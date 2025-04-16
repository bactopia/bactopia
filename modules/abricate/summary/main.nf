// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'abricate')
options.btype = "comparative"
conda_tools   = "bioconda::abricate=1.0.1"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process ABRICATE_SUMMARY {
    tag "$meta.id"
    label 'process_low'

    ext conda_tools: "bioconda::abricate=1.0.1",
        conda_name: ext.conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-"),
        conda_env: file("${params.condadir}/${ext.conda_name}").exists() ? "${params.condadir}/${ext.conda_name}" : ext.conda_tools

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/abricate:1.0.1--ha8f3691_1' :
        'quay.io/biocontainers/abricate:1.0.1--ha8f3691_1' }"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c278c7398beb73294d78639a864352abef2931ce:ba3e6d2157eac2d38d22e62ec87675e12adb1010-0':
        'biocontainers/mulled-v2-c278c7398beb73294d78639a864352abef2931ce:ba3e6d2157eac2d38d22e62ec87675e12adb1010-0' }"

    input:
    tuple val(meta), path(reports)

    output:
    tuple val(meta), path("*.tsv"), emit: report
    path "*.{log,err}"            , emit: logs, optional: true
    path ".command.*"             , emit: nf_logs
    path "versions.yml"           , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    abricate \\
        --summary \\
        $reports > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
    END_VERSIONS
    """
}
