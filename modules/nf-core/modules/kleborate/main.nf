// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'kleborate')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process KLEBORATE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? "bioconda::kleborate=2.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kleborate:2.1.0--pyhdfd78af_1' :
        'quay.io/biocontainers/kleborate:2.1.0--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(fastas)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml",emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    kleborate \\
        $options.args \\
        --outfile ${prefix}.results.txt \\
        --assemblies $fastas

    cat <<-END_VERSIONS > versions.yml
    kleborate:
        kleborate: \$( echo \$(kleborate --version | sed 's/Kleborate v//;'))
    END_VERSIONS
    """
}
