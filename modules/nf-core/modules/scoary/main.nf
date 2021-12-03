// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'scoary')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process SCOARY {
    tag "$meta.id"
    label 'process_low'
    publishDir "${publish_dir}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? "bioconda::scoary=1.6.16" : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scoary:1.6.16--py_2' :
        'quay.io/biocontainers/scoary:1.6.16--py_2' }"

    input:
    tuple val(meta), path(genes), path(traits)
    path(tree)

    output:
    tuple val(meta), path("*.csv")          , emit: csv
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs, optional: true
    path ".command.*"                       , emit: nf_logs
    path "versions.yml"                     , emit: versions

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
    scoary:
        scoary: \$( scoary --version 2>&1 )
    END_VERSIONS
    """
}
