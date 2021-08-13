include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PROKKA {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::prokka=1.14.6" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/prokka:1.14.6--pl526_0"
    } else {
        container "quay.io/biocontainers/prokka:1.14.6--pl526_0"
    }

    input:
    tuple val(meta), path(fasta)
    path proteins
    path prodigal_tf

    output:
    tuple val(meta), path("${prefix}/*.gff"), emit: gff
    tuple val(meta), path("${prefix}/*.gbk"), emit: gbk
    tuple val(meta), path("${prefix}/*.fna"), emit: fna
    tuple val(meta), path("${prefix}/*.faa"), emit: faa
    tuple val(meta), path("${prefix}/*.ffn"), emit: ffn
    tuple val(meta), path("${prefix}/*.sqn"), emit: sqn
    tuple val(meta), path("${prefix}/*.fsa"), emit: fsa
    tuple val(meta), path("${prefix}/*.tbl"), emit: tbl
    tuple val(meta), path("${prefix}/*.err"), emit: err
    tuple val(meta), path("${prefix}/*.log"), emit: log
    tuple val(meta), path("${prefix}/*.txt"), emit: txt
    tuple val(meta), path("${prefix}/*.tsv"), emit: tsv
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    def prodigal_opt = prodigal_tf ? "--prodigaltf ${prodigal_tf[0]}" : ""
    """
    prokka \\
        $options.args \\
        --cpus $task.cpus \\
        --prefix $prefix \\
        $proteins_opt \\
        $prodigal_tf \\
        $fasta

    echo \$(prokka --version 2>&1) | sed 's/^.*prokka //' > ${software}.version.txt
    """
}
