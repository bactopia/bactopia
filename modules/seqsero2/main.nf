process SEQSERO2 {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}" :
        "${task.ext.docker}" }"

    input:
    tuple val(meta), path(seqs)

    output:
    tuple val(meta), path("results/*_log.txt")   , emit: log
    tuple val(meta), path("results/*_result.tsv"), emit: tsv
    tuple val(meta), path("results/*_result.txt"), emit: txt
    path "versions.yml", emit: versions
    path ".command.begin", emit: begin
    path ".command.err", emit: err
    path ".command.log", emit: log_file
    path ".command.out", emit: out
    path ".command.run", emit: run
    path ".command.sh", emit: sh
    path ".command.trace", emit: trace

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
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
