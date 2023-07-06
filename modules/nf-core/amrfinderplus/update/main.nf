// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.containsKey("options") ? params.options : [:], 'amrfinderplus_update')
conda_tools = "bioconda::ncbi-amrfinderplus=3.11.14"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process AMRFINDERPLUS_UPDATE {
    tag "update"
    label 'process_low'
    storeDir params.datasets_cache
    publishDir params.datasets_cache

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus:3.11.14--h283d18e_1' :
        'quay.io/biocontainers/ncbi-amrfinderplus:3.11.14--h283d18e_1' }"

    output:
    path "amrfinderplus.tar.gz", emit: db

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
