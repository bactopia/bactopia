// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'

params.options = [:]
options        = initOptions(params.options, 'ectyper')
publish_dir    = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process ECTYPER {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${publish_dir}/${meta.id}",  mode: params.publish_dir_mode,  overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? "bioconda::ectyper=1.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ectyper:1.0.0--pyhdfd78af_1' :
        'quay.io/biocontainers/ectyper:1.0.0--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val(meta), path("*.txt"), emit: txt
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    ectyper \\
        $options.args \\
        --cores $task.cpus \\
        --output ./ \\
        --input $fasta_name
    mv output.tsv ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    ectyper:
        ectyper: \$(echo \$(ectyper --version 2>&1)  | sed 's/.*ectyper //; s/ .*\$//')
    END_VERSIONS
    """
}
