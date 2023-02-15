// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.containsKey("options") ? params.options : [:], 'amrfinderplus_update')
conda_tools = "bioconda::ncbi-amrfinderplus=3.10.45"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process AMRFINDERPLUS_UPDATE {
    tag "update"
    label 'process_low'
    publishDir params.outdir, mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, prefix:prefix, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus%3A3.10.45--h6e70893_0' :
        'quay.io/biocontainers/ncbi-amrfinderplus:3.10.45--h6e70893_0' }"

    output:
    path "amrfinderplus.tar.gz", emit: db
    path "*.{log,err}"       , emit: logs, optional: true
    path ".command.*"        , emit: nf_logs
    path "versions.yml"      , emit: versions

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
