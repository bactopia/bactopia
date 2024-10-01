// Import generic module functions
include { initOptions; saveFiles } from '../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'pasty')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::pasty=2.2.1" 
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process PASTY {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pasty:2.2.1--hdfd78af_0' :
        'quay.io/biocontainers/pasty:2.2.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.tsv")        , emit: tsv
    tuple val(meta), path("${prefix}.blastn.tsv") , emit: blast
    tuple val(meta), path("${prefix}.details.tsv"), emit: details
    path "*.{log,err}" , emit: logs, optional: true
    path ".command.*"  , emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    pasty \\
        $options.args \\
        --prefix $prefix \\
        --input $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pasty: \$(echo \$(pasty --version 2>&1) | sed 's/.*pasty, version //;s/ .*\$//' )
        camlhmp: \$(echo \$(pasty --version 2>&1) | sed 's/.*camlhmp, version //;s/ schema.*\$//' )
    END_VERSIONS
    """
}
