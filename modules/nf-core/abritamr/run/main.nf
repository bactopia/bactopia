// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'abritamr')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::abritamr=1.0.19"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process ABRITAMR_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/abritamr:1.0.19--pyhdfd78af_0' :
        'quay.io/biocontainers/abritamr:1.0.19--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.summary_matches.txt")  , emit: matches
    tuple val(meta), path("${prefix}.summary_partials.txt") , emit: partials
    tuple val(meta), path("${prefix}.summary_virulence.txt"), emit: virulence
    tuple val(meta), path("${prefix}.amrfinder.out")        , emit: amrfinder
    tuple val(meta), path("${prefix}.abritamr.txt")         , emit: summary, optional: true
    path "*.{log,err}"                                      , emit: logs, optional: true
    path ".command.*"                                       , emit: nf_logs
    path "versions.yml"                                     , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    abritamr run \\
        --contigs $fasta_name \\
        --prefix results \\
        $options.args \\
        --jobs $task.cpus

    # Rename output files to prevent name collisions
    mv results/summary_matches.txt ./${prefix}.summary_matches.txt
    mv results/summary_partials.txt ./${prefix}.summary_partials.txt
    mv results/summary_virulence.txt ./${prefix}.summary_virulence.txt
    mv results/amrfinder.out ./${prefix}.amrfinder.out
    if [ -f results/abritamr.txt ]; then
        # This file is not always present
        mv results/abritamr.txt ./${prefix}.abritamr.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abritamr: \$(echo \$(abritamr --version 2>&1) | sed 's/^.*abritamr //' ))
    END_VERSIONS
    """
}
