// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.containsKey("options") ? params.options : [:], 'sra-human-scrubber')
conda_tools = "bioconda::sra-human-scrubber=2.1.0"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

VERSION = '2.1.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
process SRAHUMANSCRUBBER_INITDB {
    label 'process_single'
    storeDir params.datasets_cache
    publishDir params.datasets_cache
    
    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-human-scrubber:2.1.0--hdfd78af_0' :
        'quay.io/biocontainers/sra-human-scrubber:2.1.0--hdfd78af_0' }"

    output:
    path "*.human_filter.db", emit: db
    path "*.{log,err}"      , emit: logs, optional: true
    path ".command.*"       , emit: nf_logs
    path "versions.yml"     , emit: versions

    script:
    prefix = "sra-human-scrubber-initdb"
    """
    DBVERSION=\$(curl "https://ftp.ncbi.nlm.nih.gov/sra/dbs/human_filter/current/version.txt")
    curl -f "https://ftp.ncbi.nlm.nih.gov/sra/dbs/human_filter/\${DBVERSION}.human_filter.db" -o "\${DBVERSION}.human_filter.db"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sra-human-scrubber: $VERSION
        sra-human-scrubber-db: \$DBVERSION
    END_VERSIONS
    """
}
