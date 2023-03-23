// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.containsKey("options") ? params.options : [:], 'bakta')
conda_tools = "bioconda::bakta=1.6.0"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process BAKTA_DOWNLOAD {
    label 'process_low'
    label 'process_long'
    publishDir "${params.bakta_db}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }
    
    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bakta:1.7.0--pyhdfd78af_1' :
        'quay.io/biocontainers/bakta:1.7.0--pyhdfd78af_1' }"

    output:
    path "bakta/*"     , emit: db, optional: true
    path "bakta.tar.gz", emit: db_tarball, optional: true
    path "*.{log,err}" , emit: logs, optional: true
    path ".command.*"  , emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = "bakta"
    """
    bakta_db \\
        download \\
        --output bakta

    if [ "!{params.bakta_save_as_tarball}" == "true" ]; then
        tar -czf bakta.tar.gz bakta/
        rm -rf bakta/
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$( echo \$(bakta --version 2>&1) | sed 's/^.*bakta //' )
    END_VERSIONS
    """
}
