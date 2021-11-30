// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'

options     = initOptions(params.options ? params.options : [:], 'spatyper')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process SPATYPER {
    tag "$meta.id"
    label 'process_low'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? "bioconda::spatyper=0.3.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spatyper%3A0.3.3--pyhdfd78af_3' :
        'quay.io/biocontainers/spatyper:0.3.3--pyhdfd78af_3' }"

    input:
    tuple val(meta), path(fasta)
    path repeats
    path repeat_order

    output:
    tuple val(meta), path("*.tsv")          , emit: tsv
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs, optional: true
    path ".command.*"                       , emit: nf_logs
    path "versions.yml"                     , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def input_args = repeats && repeat_order ? "-r ${repeats} -o ${repeat_order}" : ""
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    spaTyper \\
        $options.args \\
        $input_args \\
        --fasta $fasta_name \\
        --output ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    spatyper:
        spatyper: \$( echo \$(spaTyper --version 2>&1) | sed 's/^.*spaTyper //' )
    END_VERSIONS
    """
}
