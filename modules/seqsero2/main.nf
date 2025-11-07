process SEQSERO2 {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(seqs)

    output:
    tuple val(meta), path("${prefix}.tsv") , emit: tsv
    tuple val(meta), path("${prefix}.txt") , emit: txt
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
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    def is_compressed_fna = seqs[0].getName().endsWith("fna.gz") ? true : false
    def seq_name = is_compressed_fna ? seqs[0].getName().replace(".gz", "") : "${seqs}"
    """
    if [ "$is_compressed_fna" == "true" ]; then
        gzip -c -d ${seqs[0]} > $seq_name
    fi

    SeqSero2_package.py \\
        ${task.ext.args} \\
        -d supplemental/ \\
        -n $prefix \\
        -p $task.cpus \\
        -t 4 \\
        -i $seq_name

    mv supplemental/SeqSero_log.txt ./${prefix}.log
    mv supplemental/SeqSero_result.tsv ./${prefix}.tsv
    mv supplemental/SeqSero_result.txt ./${prefix}.txt

    # Cleanup
    rm -rf supplemental/ ${seq_name} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqsero2: \$( echo \$( SeqSero2_package.py --version 2>&1) | sed 's/^.*SeqSero2_package.py //' )
    END_VERSIONS
    """
}
