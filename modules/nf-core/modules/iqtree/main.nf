// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'

params.options = [:]
options        = initOptions(params.options, iqtree)
publish_dir    = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process IQTREE {
    tag "$alignment"
    label 'process_medium'
    publishDir "${publish_dir}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? 'bioconda::iqtree=2.1.4_beta' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/iqtree:2.1.4_beta--hdcc8f71_0' :
        'quay.io/biocontainers/iqtree:2.1.4_beta--hdcc8f71_0' }"

    input:
    path alignment
    val constant_sites

    output:
    path "*.treefile", emit: phylogeny
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    script:
    def fconst_args = constant_sites ? "-fconst $constant_sites" : ''
    def memory      = task.memory.toString().replaceAll(' ', '')
    """
    iqtree \\
        $fconst_args \\
        $options.args \\
        -s $alignment \\
        -nt AUTO \\
        -ntmax $task.cpus \\
        -mem $memory \\

    cat <<-END_VERSIONS > versions.yml
    iqtree:
        iqtree: \$(echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
    END_VERSIONS
    """
}
