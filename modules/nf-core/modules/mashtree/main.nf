// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../../../../lib/nf/functions'

params.options = [:]
options        = initOptions(params.options)
publish_dir    = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process MASHTREE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${publish_dir}",
        mode: params.publish_dir_mode,
        overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, process_name:getSoftwareName(task.process, options.full_software_name), is_module: options.is_module) }

    conda (params.enable_conda ? "bioconda::mashtree=1.2.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mashtree:1.2.0--pl526h516909a_0"
    } else {
        container "quay.io/biocontainers/mashtree:1.2.0--pl526h516909a_0"
    }

    input:
    tuple val(meta), path(seqs)

    output:
    tuple val(meta), path("*.dnd"), emit: tree
    tuple val(meta), path("*.tsv"), emit: matrix
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml",emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    mashtree \\
        $options.args \\
        --numcpus $task.cpus \\
        --outmatrix ${prefix}.tsv \\
        --outtree ${prefix}.dnd \\
        $seqs

    cat <<-END_VERSIONS > versions.yml
    mashtree:
        mash: \$(echo \$(mash --version 2>&1))
        mashtree: \$( echo \$( mashtree --version 2>&1 ) | sed 's/^.*Mashtree //' )
    END_VERSIONS
    """
}
