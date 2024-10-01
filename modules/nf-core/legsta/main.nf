// Import generic module functions
include { initOptions; saveFiles } from '../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'legsta')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::legsta=0.5.1"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process LEGSTA {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/legsta:0.5.1--hdfd78af_2' :
        'quay.io/biocontainers/legsta:0.5.1--hdfd78af_2' }"

    input:
    tuple val(meta), path(seqs)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "*.{log,err}"            , emit: logs, optional: true
    path ".command.*"             , emit: nf_logs
    path "versions.yml"           , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    legsta \\
        $options.args \\
        $seqs | sed 's/.fna//; s/.gz//' > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        legsta: \$(echo \$(legsta --version 2>&1) | sed 's/^.*legsta //; s/ .*\$//;')
    END_VERSIONS
    """
}
