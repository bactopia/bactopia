// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ISMAPPER {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::ismapper=2.0.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ismapper:2.0.2--pyhdfd78af_1"
    } else {
        container "quay.io/biocontainers/ismapper:2.0.2--pyhdfd78af_1"
    }

    input:
    tuple val(meta), path(reads), path(reference), path(query)

    output:
    tuple val(meta), path("results/*"), emit: results
    path "versions.yml"               , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    ismap \\
        $options.args \\
        --t $task.cpus \\
        --output_dir results \\
        --queries $query \\
        --reference $reference \\
        --reads $reads

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$( ismap --version 2>&1 ) | sed 's/^.*ismap //' )
    END_VERSIONS
    """
}
