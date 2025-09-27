process SYLPH_PROFILE {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(reads)
    path db

    output:
    tuple val(meta), path("${prefix}.tsv") , emit: tsv
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
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    def query_reads = meta.single_end ? "${reads}" : "--first-pairs ${reads[0]} --second-pairs ${reads[1]}"
    """
    sylph \\
        profile \\
        $db \\
        $query_reads \\
        -t ${task.cpus} \\
        ${task.ext.args} \\
        --output-file ${prefix}.original.tsv

    # Remove the "fasta.gz" from sample names in output
    if [ "${meta.single_end}" == "true" ]; then
        sed 's/^${prefix}.fastq.gz/${prefix}/' ${prefix}.original.tsv > ${prefix}.tsv
    else
        sed 's/^${prefix}_R1.fastq.gz/${prefix}/' ${prefix}.original.tsv > ${prefix}.tsv
    fi
    rm ${prefix}.original.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(echo \$(sylph --version 2>&1) | sed 's/^.*sylph //;s/ .*\$//')
    END_VERSIONS
    """
}
