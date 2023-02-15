// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../lib/nf/functions'
RESOURCES     = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options       = initOptions(params.containsKey("options") ? params.options : [:], 'iqtree')
options.btype = options.btype ?: "comparative"
conda_tools   = "bioconda::iqtree=2.2.0.3"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process IQTREE {
    tag "$prefix"
    label 'process_medium'
    label 'process_long'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/iqtree:2.2.0.3--hb97b32f_1' :
        'quay.io/biocontainers/iqtree:2.2.0.3--hb97b32f_1' }"

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("${prefix}*")                         , emit: results
    tuple val(meta), path("${prefix}.treefile")                 , emit: phylogeny
    tuple val(meta), path(alignment), path("${prefix}.treefile"), emit: aln_tree
    path "*.{log,err}"                                          , emit: logs, optional: true
    path ".command.*"                                           , emit: nf_logs
    path "versions.yml"                                         , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def memory = task.memory.toString().replaceAll(' ', '')
    """
    iqtree \\
        $options.args \\
        -s $alignment \\
        -nt $task.cpus \\
        -ntmax $task.cpus \\
        -pre $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
    END_VERSIONS
    """
}
