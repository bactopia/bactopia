// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../lib/nf/functions'
RESOURCES     = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options       = initOptions(params.containsKey("options") ? params.options : [:], 'agrvate')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::agrvate=1.0.2"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process AGRVATE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agrvate:1.0.2--hdfd78af_0' :
        'quay.io/biocontainers/agrvate:1.0.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("results/${meta.id}-summary.tab"), emit: summary
    tuple val(meta), path("results/*")                     , emit: results_dir
    path "*.{log,err}"                                     , emit: logs, optional: true
    path ".command.*"                                      , emit: nf_logs
    path "versions.yml"                                    , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    agrvate \\
        $options.args \\
        -i $fasta_name

    mv ${meta.id}-results/ results/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agrvate: \$(echo \$(agrvate -v 2>&1) | sed 's/agrvate v//;')
    END_VERSIONS
    """
}
