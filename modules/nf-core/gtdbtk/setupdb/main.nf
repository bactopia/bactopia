// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions([:], 'gtdb')
conda_tools = "bioconda::gtdbtk=2.1.1"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process GTDBTK_SETUPDB {
    label 'process_high'
    publishDir "${params.gtdb}", mode: params.publish_dir_mode, overwrite: true,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gtdbtk:2.1.1--pyhdfd78af_1' :
        'quay.io/biocontainers/gtdbtk:2.1.1--pyhdfd78af_1' }"

    output:
    path("gtdbtk/*")     , emit: db, optional: true
    path("gtdbtk.tar.gz"), emit: db, optional: true
    path "*.{log,err}"   , emit: logs, optional: true
    path ".command.*"    , emit: nf_logs
    path "versions.yml"  , emit: versions

    script:
    """
    export GTDBTK_DATA_PATH="./gtdbtk"
    mkdir ./gtdbtk
    download-db.sh ./gtdbtk
    gtdbtk check_install && touch gtdb-setup.txt

    if [ "!{params.gtdb_save_as_tarball}" == "true" ]; then
        tar -czf gtdbtk.tar.gz gtdbtk/
        rm -rf gtdbtk/
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo \$(gtdbtk --version -v 2>&1) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """
}
