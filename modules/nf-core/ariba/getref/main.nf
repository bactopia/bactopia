// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.containsKey("options") ? params.options : [:], 'ariba')
conda_tools = "bioconda::ariba=2.14.6 bcftools=1.14 pysam=0.18.0"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process ARIBA_GETREF {
    tag "$db_name"
    label 'process_low'
    storeDir params.datasets_cache
    publishDir params.datasets_cache

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ariba:2.14.6--py39h67e14b5_3' :
        'quay.io/biocontainers/ariba:2.14.6--py39h67e14b5_3' }"

    input:
    val(db_name)

    output:
    path("ariba-${db_name}.tar.gz"), emit: db

    script:
    """
    # Download, format database, and tarball it
    ariba \\
        getref \\
        ${db_name} \\
        ${db_name}

    ariba \\
        prepareref \\
        -f ${db_name}.fa \\
        -m ${db_name}.tsv \\
        ${db_name}

    mv ${db_name}/ ariba-${db_name}/
    tar -zcvf ariba-${db_name}.tar.gz ariba-${db_name}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ariba:  \$(echo \$(ariba version 2>&1) | sed 's/^.*ARIBA version: //;s/ .*\$//')
    END_VERSIONS
    """
}
