// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.containsKey("options") ? params.options : [:], 'eggnog_download')
conda_tools = "bioconda::eggnog-mapper=2.1.11"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process EGGNOG_DOWNLOAD {
    label 'process_low'
    label 'process_long'
    publishDir "${params.eggnog_db}", mode: params.publish_dir_mode, overwrite: true

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.11--pyhdfd78af_0' :
        'quay.io/biocontainers/eggnog-mapper:2.1.11--pyhdfd78af_0' }"

    output:
    path("eggnog/")     , emit: db, optional: true
    path("eggnog.tar.gz"), emit: db_tarball, optional: true
    path "*.{log,err}"   , emit: logs, optional: true
    path ".command.*"    , emit: nf_logs
    path "versions.yml"  , emit: versions

    script:
    """
    mkdir eggnog
    download_eggnog_data.py \\
        $options.args \\
        -y \\
        --data_dir eggnog/

    if [ "${params.eggnog_save_as_tarball}" == "true" ]; then
        tar -czf eggnog.tar.gz eggnog/
        rm -rf eggnog/
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$( echo \$(emapper.py --version 2>&1)| sed 's/.* emapper-//;s/ .*//')
    END_VERSIONS
    """
}
