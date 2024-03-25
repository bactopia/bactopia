// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES     = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options       = initOptions(params.containsKey("options") ? params.options : [:], 'custom_wget')
options.btype = options.btype ?: "comparative"
conda_tools   = "bioconda::gnu-wget=1.18"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process CUSTOM_WGET {
    label 'process_low'
    storeDir params.datasets_cache
    publishDir params.datasets_cache

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9' :
        'quay.io/biocontainers/gnu-wget:1.18--h36e9172_9' }"

    output:
    path "${options.args2}", emit: download
    path "*.{log,err}"   , emit: logs, optional: true
    path ".command.*"    , emit: nf_logs
    path "versions.yml"  ,emit: versions

    script:
    """
    wget $options.args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo \$( wget --version 2>&1) | sed "s/GNU Wget //;s/ .*//;" )
    END_VERSIONS
    """
}
