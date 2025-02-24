// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'sylph')
options.btype = "tools"
conda_tools   = "bioconda::sylph=0.8.0"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process SYLPH_PROFILE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sylph:0.8.0--h919a2d8_0' :
        'quay.io/biocontainers/sylph:0.8.0--h919a2d8_0' }"

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def query_reads = meta.single_end ? "${reads}" : "--first-pairs ${reads[0]} --second-pairs ${reads[1]}"
    """
    sylph \\
        profile \\
        $reference \\
        $query_reads \\
        -t ${task.cpus} \\
        $options.args \\
        --output-file ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(echo \$(sylph --version 2>&1) | sed 's/^.*sylph //;s/ .*\$//')
    END_VERSIONS
    """
}
