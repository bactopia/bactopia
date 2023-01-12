// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.containsKey("options") ? params.options : [:], 'amrfinderplus_update')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::ncbi-amrfinderplus=3.10.45"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process AMRFINDERPLUS_UPDATE {
    tag "update"
    label 'process_low'
    publishDir "${publish_dir}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus%3A3.10.45--h6e70893_0' :
        'quay.io/biocontainers/ncbi-amrfinderplus:3.10.45--h6e70893_0' }"

    output:
    path "amrfinderdb.tar.gz", emit: db
    path "*.{log,err}"       , emit: logs, optional: true
    path ".command.*"        , emit: nf_logs
    path "versions.yml"      , emit: versions

    script:
    """
    mkdir amrfinderdb
    amrfinder_update -d amrfinderdb
    tar czvf amrfinderdb.tar.gz -C amrfinderdb/\$(readlink amrfinderdb/latest) ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
    END_VERSIONS
    """
}
