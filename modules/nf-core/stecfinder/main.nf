// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../lib/nf/functions'
RESOURCES     = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options       = initOptions(params.containsKey("options") ? params.options : [:], 'stecfinder')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::stecfinder=1.1.0"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process STECFINDER {
    tag "$meta.id"
    label 'process_low'
    publishDir params.outdir, mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, prefix:prefix, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stecfinder:1.1.0--pyhdfd78af_0':
        'quay.io/biocontainers/stecfinder:1.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(seqs)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "*.{log,err}"            , emit: logs, optional: true
    path ".command.*"             , emit: nf_logs
    path "versions.yml"           , emit: versions

    script:
    def prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def is_compressed = meta.is_compressed && !params.stecfinder_use_reads ? true : false
    def seq_name = is_compressed ? seqs[0].getName().replace(".gz", "") : seqs
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d $seqs > $seq_name
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
