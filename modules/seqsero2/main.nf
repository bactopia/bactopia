process SEQSERO2 {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(seqs)

    output:
    tuple val(meta), path("results/*_log.txt")   , emit: log
    tuple val(meta), path("results/*_result.tsv"), emit: tsv
    tuple val(meta), path("results/*_result.txt"), emit: txt
    tuple val(meta), path("*.{log,err}")   , emit: logs, optional: true
    tuple val(meta), path(".command.begin"), emit: nf_begin
    tuple val(meta), path(".command.err")  , emit: nf_err
    tuple val(meta), path(".command.log")  , emit: nf_log_file
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
    def is_compressed_fna = seqs[0].getName().endsWith("fna.gz") ? true : false
    def seq_name = is_compressed_fna ? seqs[0].getName().replace(".gz", "") : "${seqs}"
    """
    if [ "$is_compressed_fna" == "true" ]; then
        gzip -c -d ${seqs[0]} > $seq_name
    fi

    SeqSero2_package.py \\
        $args \\
        -d results/ \\
        -n $prefix \\
        -p $task.cpus \\
        -i $seq_name

    mv results/SeqSero_log.txt results/${prefix}_log.txt
    mv results/SeqSero_result.tsv results/${prefix}_result.tsv
    mv results/SeqSero_result.txt results/${prefix}_result.txt

    # Cleanup
    rm -rf results/${seq_name} ${seq_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqsero2: \$( echo \$( SeqSero2_package.py --version 2>&1) | sed 's/^.*SeqSero2_package.py //' )
    END_VERSIONS
    """
}
