// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../lib/nf/functions'
RESOURCES     = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options       = initOptions(params.containsKey("options") ? params.options : [:], 'gamma')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::gamma=2.2"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools
def VERSION   = '2.1' // Version information not provided by tool on CLI

process GAMMA {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gamma:2.2--hdfd78af_0' :
        'quay.io/biocontainers/gamma:2.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    tuple val(meta), path("*.gamma")                , emit: gamma
    tuple val(meta), path("*.psl")                  , emit: psl
    tuple val(meta), path("*.gff")  , optional:true , emit: gff
    tuple val(meta), path("*.fasta"), optional:true , emit: fasta
    path "*.{log,err}"            , emit: logs, optional: true
    path ".command.*"             , emit: nf_logs
    path "versions.yml"           , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    GAMMA.py \\
        $options.args \\
        $fasta_name \\
        $db \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gamma: $VERSION
    END_VERSIONS
    """
}
