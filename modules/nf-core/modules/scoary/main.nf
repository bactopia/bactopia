// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SCOARY {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::scoary=1.6.16" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/scoary:1.6.16--py_2"
    } else {
        container "quay.io/biocontainers/scoary:1.6.16--py_2"
    }

    input:
    tuple val(meta), path(genes), path(traits)
    path(tree)

    output:
    tuple val(meta), path("*.csv"), emit: csv
    path "versions.yml"           , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def newick_tree = tree ? "-n ${tree}" : ""
    """
    scoary \\
        $options.args \\
        --no-time \\
        --threads $task.cpus \\
        --traits $traits \\
        --genes $genes

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( scoary --version 2>&1 )
    END_VERSIONS
    """
}
