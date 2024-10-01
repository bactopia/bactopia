// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'rgi')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::rgi=6.0.3"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process RGI_MAIN {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rgi:6.0.3--pyha8f3691_0' :
        'quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.json"), emit: json, optional: true
    tuple val(meta), path("*.txt") , emit: tsv
    path "*.{log,err}"             , emit: logs, optional: true
    path ".command.*"              , emit: nf_logs
    path "versions.yml"            , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    rgi \\
        main \\
        $options.args \\
        --clean \\
        --data wgs \\
        --num_threads $task.cpus \\
        --output_file $prefix \\
        --input_sequence $fasta

    # Remove empty json files
    if grep "^{}\$" ${prefix}.json; then
        rm ${prefix}.json
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(rgi main --version)
    END_VERSIONS
    """
}
