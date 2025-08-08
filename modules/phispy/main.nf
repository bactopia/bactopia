process PHISPY {
    tag "$meta.id"
    label 'process_medium'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}" :
        "${task.ext.docker}" }"

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
    path "versions.yml"                                                 , emit: versions
    path ".command.begin"                                               , emit: begin
    path ".command.err"                                                 , emit: err
    path ".command.log"                                                 , emit: log
    path ".command.out"                                                 , emit: out
    path ".command.run"                                                 , emit: run
    path ".command.sh"                                                  , emit: sh
    path ".command.trace"                                               , emit: trace

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir results/
    PhiSpy.py \\
        $args \\
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
