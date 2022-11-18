// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'rgi')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::rgi=6.0.1"
conda_env   = file("${params.condadir}/rgi").exists() ? "${params.condadir}/rgi" : conda_tools

process RGI_MAIN {
    tag "$meta.id"
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rgi:6.0.1--pyha8f3691_1' :
        'quay.io/biocontainers/rgi:6.0.1--pyha8f3691_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.txt") , emit: tsv
    path "*.{log,err}"             , emit: logs, optional: true
    path ".command.*"              , emit: nf_logs
    path "versions.yml"            , emit: versions

    script:
    def prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    rgi \\
        main \\
        $options.args \\
        --clean \\
        --data wgs \\
        --num_threads $task.cpus \\
        --output_file $prefix \\
        --input_sequence $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(rgi main --version)
    END_VERSIONS
    """
}
