// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process IQTREE {
    tag "$alignment"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? 'bioconda::iqtree=2.1.4_beta' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/iqtree:2.1.4_beta--hdcc8f71_0"
    } else {
        container "quay.io/biocontainers/iqtree:2.1.4_beta--hdcc8f71_0"
    }

    input:
    path alignment
    val constant_sites

    output:
    path "*.treefile",    emit: phylogeny
    path "*.version.txt", emit: version

    script:
    def software    = getSoftwareName(task.process)
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

    echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version \\([0-9\\.]*\\) .*\$/\\1/' > ${software}.version.txt
    """
}
