// Import generic module functions
include { initOptions; saveFiles } from '../../../../../lib/nf/functions'

options     = initOptions(params.options ? params.options : [:], 'csvtk_concat')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir

process CSVTK_CONCAT {
    tag "$meta.id"
    label 'process_low'
    publishDir "${publish_dir}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? "bioconda::csvtk=0.23.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/csvtk:0.23.0--h9ee0642_0' :
        'quay.io/biocontainers/csvtk:0.23.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(csv)
    val in_format
    val out_format

    output:
    tuple val(meta), path("${prefix}.${out_extension}"), emit: csv
    path "*.{stdout.txt,stderr.txt,log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml",emit: versions

    script:
    prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def delimiter = in_format == "tsv" ? "\t" : (in_format == "csv" ? "," : in_format)
    def out_delimiter = out_format == "tsv" ? "\t" : (out_format == "csv" ? "," : out_format)
    out_extension = out_format == "tsv" ? 'tsv' : 'csv'
    """
    csvtk \\
        concat \\
        $options.args \\
        --num-cpus $task.cpus \\
        --delimiter "${delimiter}" \\
        --out-delimiter "${out_delimiter}" \\
        --out-file ${prefix}.${out_extension} \\
        $csv

    cat <<-END_VERSIONS > versions.yml
    csvtk_concat:
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """
}
