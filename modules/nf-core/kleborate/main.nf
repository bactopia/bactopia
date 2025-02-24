// Import generic module functions
include { initOptions; saveFiles } from '../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'kleborate')
options.btype = "tools"
conda_tools   = "bioconda::kleborate=3.1.3"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process KLEBORATE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kleborate:3.1.3--pyhdfd78af_0' :
        'quay.io/biocontainers/kleborate:3.1.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastas)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml",emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    mkdir results/
    kleborate \\
        $options.args \\
        --outdir results/ \\
        --assemblies $fastas

    # Rename output file to include the prefix name
    find results/ -name "*output.txt" -print0 | while read -d \$'\0' file; do mv "\$file" "${prefix}.txt"; done

    # cleanup
    rm -rf results/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kleborate: \$( echo \$(kleborate --version | sed 's/Kleborate v//;'))
    END_VERSIONS
    """
}
