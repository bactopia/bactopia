// Import generic module functions
include { initOptions; saveFiles } from '../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'stecfinder')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::stecfinder=1.1.1"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process STECFINDER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stecfinder:1.1.1--pyhdfd78af_0' :
        'quay.io/biocontainers/stecfinder:1.1.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(reads)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "*.{log,err}"            , emit: logs, optional: true
    path ".command.*"             , emit: nf_logs
    path "versions.yml"           , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_compressed = meta.is_compressed && !params.stecfinder_use_reads ? true : false
    def seq_name = is_compressed ? fasta.getName().replace(".gz", "") : reads
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d $fasta > $seq_name
    fi

    stecfinder \\
        -i $seq_name \\
        $options.args \\
        -t $task.cpus > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stecfinder: \$(echo \$(stecfinder --version 2>&1) | sed 's/^.*STECFinder version: //;' )
    END_VERSIONS
    """
}
