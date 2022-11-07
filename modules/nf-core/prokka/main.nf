// Import generic module functions
include { get_resources; initOptions; saveFiles } from '../../../lib/nf/functions'
RESOURCES   = get_resources(workflow.profile, params.max_memory, params.max_cpus)
options     = initOptions(params.options ? params.options : [:], 'prokka')
publish_dir = params.is_subworkflow ? "${params.outdir}/bactopia-tools/${params.wf}/${params.run_name}" : params.outdir
conda_tools = "bioconda::prokka=1.14.6 perl-bioperl=1.7.2"
conda_name  = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env   = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process PROKKA {
    tag "$meta.id"
    label 'process_low'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/prokka:1.14.6--pl526_0' :
        'quay.io/biocontainers/prokka:1.14.6--pl526_0' }"

    input:
    tuple val(meta), path(fasta)
    path proteins
    path prodigal_tf

    output:
    tuple val(meta), path("${prefix}/*.{ffn,ffn.gz}"), path("${prefix}/*.{faa,faa.gz}"), emit: annotations
    tuple val(meta), path("${prefix}/*.{gff,gff.gz}"), emit: gff
    tuple val(meta), path("${prefix}/*.{gbk,gbk.gz}"), emit: gbk
    tuple val(meta), path("${prefix}/*.{fna,fna.gz}"), emit: fna
    tuple val(meta), path("${prefix}/*.{faa,faa.gz}"), emit: faa
    tuple val(meta), path("${prefix}/*.{ffn,ffn.gz}"), emit: ffn
    tuple val(meta), path("${prefix}/*.{sqn,sqn.gz}"), emit: sqn
    tuple val(meta), path("${prefix}/*.{fsa,fsa.gz}"), emit: fsa
    tuple val(meta), path("${prefix}/*.{tbl,tbl.gz}"), emit: tbl
    tuple val(meta), path("${prefix}/*.txt"), emit: txt
    tuple val(meta), path("${prefix}/*.tsv"), emit: tsv
    path "*.{log,err}", emit: logs, optional: true
    path ".command.*", emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    def prodigal_opt = prodigal_tf ? "--prodigaltf ${prodigal_tf[0]}" : ""
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    prokka \\
        $options.args \\
        --cpus $task.cpus \\
        --prefix $prefix \\
        $proteins_opt \\
        $prodigal_tf \\
        $fasta_name

    if [[ "${params.skip_compression}" == "false" ]]; then
        gzip ${prefix}/*.gff
        gzip ${prefix}/*.gbk
        gzip ${prefix}/*.fna
        gzip ${prefix}/*.faa
        gzip ${prefix}/*.ffn
        gzip ${prefix}/*.sqn
        gzip ${prefix}/*.fsa
        gzip ${prefix}/*.tbl
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prokka: \$( echo \$(prokka --version 2>&1) | sed 's/^.*prokka //')
    END_VERSIONS
    """
}
