// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options        = initOptions(params.containsKey("options") ? params.options : [:], 'defensefinder_update')
conda_tools    = "bioconda::defense-finder=1.2.2"
conda_name     = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env      = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools
DF_VERSION     = "1.3.0"
DF_MODELS_VERSION = "1.3.0"
CASFINDER_VERSION = "3.1.0"

process DEFENSEFINDER_UPDATE {
    tag "update"
    label 'process_low'
    storeDir params.datasets_cache
    publishDir params.datasets_cache

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/defense-finder:1.2.2--pyhdfd78af_0' :
        'quay.io/biocontainers/defense-finder:1.2.2--pyhdfd78af_0' }"

    output:
    path "defense-finder-models-${DF_MODELS_VERSION}.tar", emit: db

    script:
    prefix = "defense-finder"
    """
    mkdir models
    wget \\
        -O models/defense-finder-models-v${DF_MODELS_VERSION}.tar.gz \\
        https://github.com/mdmparis/defense-finder-models/archive/refs/tags/${DF_MODELS_VERSION}.tar.gz

    wget \\
        -O models/CasFinder-${CASFINDER_VERSION}.tar.gz \\
        https://github.com/macsy-models/CasFinder/archive/refs/tags/${CASFINDER_VERSION}.tar.gz

    tar -cvf defense-finder-models-${DF_MODELS_VERSION}.tar models/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        defense-finder: ${DF_VERSION}
        defense-finder-models: ${DF_MODELS_VERSION}
        casfinder-models: ${CASFINDER_VERSION}
    END_VERSIONS
    """
}
