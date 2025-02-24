// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options     = initOptions(params.containsKey("options") ? params.options : [:], 'bakta')
conda_tools = "bioconda::bakta=1.10.4"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process BAKTA_DOWNLOAD {
    label 'process_low'
    label 'process_long'

    publishDir "${params.bakta_db}", mode: params.publish_dir_mode, overwrite: true
    
    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bakta:1.10.4--pyhdfd78af_0' :
        'quay.io/biocontainers/bakta:1.10.4--pyhdfd78af_0' }"

    output:
    path "bakta-${params.bakta_db_type}/*"     , emit: db, optional: true
    path "bakta-${params.bakta_db_type}.tar.gz", emit: db_tarball, optional: true
    path "*.{log,err}" , emit: logs, optional: true
    path ".command.*"  , emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = "bakta"
    """
    bakta_db \\
        download \\
        --type "${params.bakta_db_type}" \\
        --output bakta

    if [ "${params.bakta_save_as_tarball}" == "true" ]; then
        tar -czf bakta-${params.bakta_db_type}.tar.gz bakta/
        rm -rf bakta/
    else
        mv bakta/ bakta-${params.bakta_db_type}/
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$( echo \$(bakta --version 2>&1) | sed 's/^.*bakta //' )
    END_VERSIONS
    """
}
