// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'

options     = initOptions(params.options ? params.options : [:], 'lissero')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process LISSERO {
    tag "$meta.id"
    label 'process_low'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? "bioconda::lissero=0.4.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lissero:0.4.9--py_0' :
        'quay.io/biocontainers/lissero:0.4.9--py_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv")          , emit: tsv
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs, optional: true
    path ".command.*"                       , emit: nf_logs
    path "versions.yml"                     , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    lissero \\
        $options.args \\
        $fasta_name \\
        > ${prefix}.tsv
    sed -i 's/^.*${fasta_name}/${fasta_name}/' ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    lissero:
        lissero: \$( echo \$(lissero --version 2>&1) | sed 's/^.*LisSero //' )
    END_VERSIONS
    """
}
