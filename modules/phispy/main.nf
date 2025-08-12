process PHISPY {
    tag "$meta.id"
    label 'process_medium'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

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
    tuple val(meta), path("*.{log,err}")   , emit: logs, optional: true
    tuple val(meta), path(".command.begin"), emit: nf_begin
    tuple val(meta), path(".command.err")  , emit: nf_err
    tuple val(meta), path(".command.log")  , emit: nf_log
    tuple val(meta), path(".command.out")  , emit: nf_out
    tuple val(meta), path(".command.run")  , emit: nf_run
    tuple val(meta), path(".command.sh")   , emit: nf_sh
    tuple val(meta), path(".command.trace"), emit: nf_trace
    tuple val(meta), path("versions.yml")  , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs"
    meta.process_name = task.ext.process_name
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
