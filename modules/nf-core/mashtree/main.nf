// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.containsKey("options") ? params.options : [:], 'mashtree')
options.btype = options.btype ?: "comparative"
conda_tools = "bioconda::mashtree=1.2.0"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process MASHTREE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mashtree:1.2.0--pl5321hec16e2b_1' :
        'quay.io/biocontainers/mashtree:1.2.0--pl5321hec16e2b_1' }"

    input:
    tuple val(meta), path(seqs)

    output:
    tuple val(meta), path("*.dnd"), emit: tree
    tuple val(meta), path("*.tsv"), emit: matrix
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml",emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    mashtree \\
        $options.args \\
        --numcpus $task.cpus \\
        --outmatrix ${prefix}.tsv \\
        --outtree ${prefix}.dnd \\
        $seqs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(echo \$(mash 2>&1) | sed 's/^.*Mash version //;s/ .*\$//')
        mashtree: \$( echo \$( mashtree --version 2>&1 ) | sed 's/^.*Mashtree //' )
    END_VERSIONS
    """
}
