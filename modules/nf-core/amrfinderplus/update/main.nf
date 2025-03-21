// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options     = initOptions(params.containsKey("options") ? params.options : [:], 'amrfinderplus_update')
conda_tools = "bioconda::ncbi-amrfinderplus=4.0.19"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process AMRFINDERPLUS_UPDATE {
    tag "update"
    label 'process_low'
    storeDir params.datasets_cache
    publishDir params.datasets_cache

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus:4.0.19--hf69ffd2_0' :
        'quay.io/biocontainers/ncbi-amrfinderplus:4.0.19--hf69ffd2_0' }"

    output:
    path "amrfinderplus.tar.gz", emit: db
    path "*.{log,err}" , emit: logs, optional: true
    path ".command.*"  , emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = "amrfinderplus"
    """
    mkdir amrfinderplus-temp
    amrfinder_update -d amrfinderplus-temp
    mv amrfinderplus-temp/\$(readlink amrfinderplus-temp/latest) amrfinderplus/
    tar czvf amrfinderplus.tar.gz amrfinderplus/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
        amrfinderplus-database: \$(echo \$(echo \$(amrfinder --database amrfinderplus --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev))
    END_VERSIONS
    """
}
