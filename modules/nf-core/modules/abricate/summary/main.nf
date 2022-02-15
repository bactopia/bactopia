// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'abricate')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::abricate=1.0.1"
conda_env   = file("${params.condadir}/abricate").exists() ? "${params.condadir}/abricate" : conda_tools

process ABRICATE_SUMMARY {
    tag "$meta.id"
    label 'process_low'
    publishDir "${publish_dir}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/abricate%3A1.0.1--ha8f3691_1' :
        'quay.io/biocontainers/abricate:1.0.1--ha8f3691_1' }"

    input:
    tuple val(meta), path(reports)

    output:
    tuple val(meta), path("*.txt"), emit: report
    path "*.{log,err}"            , emit: logs, optional: true
    path ".command.*"             , emit: nf_logs
    path "versions.yml"           , emit: versions

    script:
    def prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    abricate \\
        --summary \\
        $reports > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
    END_VERSIONS
    """
}
