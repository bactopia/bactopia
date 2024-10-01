// Import generic module functions
include { initOptions; saveFiles } from '../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'phispy')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::phispy=4.2.21"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

process PHISPY {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/phispy:4.2.21--py310h0dbaff4_2' :
        'quay.io/biocontainers/phispy:4.2.21--py310h0dbaff4_2' }"

    input:
    tuple val(meta), path(gbk)

    output:
    tuple val(meta), path("results/${prefix}.tsv")                     , emit: tsv
    tuple val(meta), path("results/${prefix}_prophage_information.tsv"), optional:true, emit: information
    tuple val(meta), path("results/${prefix}_bacteria.fasta")          , optional:true, emit: bacteria_fasta
    tuple val(meta), path("results/${prefix}_bacteria.gbk")            , optional:true, emit: bacteria_gbk
    tuple val(meta), path("results/${prefix}_phage.fasta")             , optional:true, emit: phage_fasta
    tuple val(meta), path("results/${prefix}_phage.gbk")               , optional:true, emit: phage_gbk
    tuple val(meta), path("results/${prefix}_prophage.gff3")           , optional:true, emit: prophage_gff
    tuple val(meta), path("results/${prefix}_prophage.tbl")            , optional:true, emit: prophage_tbl
    tuple val(meta), path("results/${prefix}_prophage.tsv")            , optional:true, emit: prophage_tsv
    path "*.{log,err}" , emit: logs, optional: true
    path ".command.*"  , emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix = options.suffix ? "${options.suffix}" : "${meta.id}"
    """
    mkdir results/
    PhiSpy.py \\
        $options.args \\
        --threads $task.cpus \\
        -p $prefix \\
        -o results \\
        $gbk

    mv results/${prefix}_prophage_coordinates.tsv results/${prefix}.tsv
    mv results/${prefix}_phispy.log ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PhiSpy: \$(echo \$(PhiSpy.py --version 2>&1))
    END_VERSIONS
    """
}
