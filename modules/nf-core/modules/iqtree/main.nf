// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'iqtree')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process IQTREE {
    tag "$prefix"
    label 'process_medium'
    publishDir "${publish_dir}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? 'bioconda::iqtree=2.1.4_beta' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/iqtree:2.1.4_beta--hdcc8f71_0' :
        'quay.io/biocontainers/iqtree:2.1.4_beta--hdcc8f71_0' }"

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("${prefix}*")                         , emit: results
    tuple val(meta), path("${prefix}.treefile")                 , emit: phylogeny
    tuple val(meta), path(alignment), path("${prefix}.treefile"), emit: aln_tree
    path "*.{stdout.txt,stderr.txt,log,err}"                    , emit: logs, optional: true
    path ".command.*"                                           , emit: nf_logs
    path "versions.yml"                                         , emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def memory = task.memory.toString().replaceAll(' ', '')
    """
    iqtree \\
        $options.args \\
        -s $alignment \\
        -nt AUTO \\
        -ntmax $task.cpus \\
        -pre $prefix

    cat <<-END_VERSIONS > versions.yml
    iqtree:
        iqtree: \$(echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
    END_VERSIONS
    """
}
