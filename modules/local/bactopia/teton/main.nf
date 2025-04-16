// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'teton-prepare')
options.btype = "tools"
conda_tools   = "bioconda::bactopia-teton=1.1.1"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process BACTOPIA_SAMPLESHEET {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bactopia-teton:1.1.1--hdfd78af_0' :
        'quay.io/biocontainers/bactopia-teton:1.1.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(classification)

    output:
    tuple val(meta), path("${prefix}.bacteria.tsv")   , emit: bacteria_tsv
    tuple val(meta), path("${prefix}.nonbacteria.tsv"), emit: nonbacteria_tsv
    tuple val(meta), path("${prefix}-sizemeup.txt")   , emit: sizemeup
    path "*.{log,err}"         , emit: logs, optional: true
    path ".command.*"          , emit: nf_logs
    path "versions.yml"        , emit: versions
    path "*-{error,merged}.txt", optional: true

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    # determine genome size and create sample sheet
    sizemeup \\
        --query $classification \\
        --prefix ${prefix}

    # create sample sheet
    teton-prepare.py \\
        ${prefix} \\
        ${prefix}-sizemeup.txt \\
        ${meta.runtype} \\
        ${meta.teton_reads} \\
        ${params.outdir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sizemeup: \$(echo \$(sizemeup --version 2>&1) | sed 's/.*sizemeup-main, version //;s/ .*\$//' )
    END_VERSIONS
    """
}
