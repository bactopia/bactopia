process SYLPH_PROFILE {
    tag "$meta.id"
    label 'process_medium'

    conda "${task.ext.conda_env}"
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sylph:0.8.0--ha6fb395_0' :
        'quay.io/biocontainers/sylph:0.8.0--ha6fb395_0' }"

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    path "*.{log,err}", emit: logs, optional: true
    path ".command.begin"   , emit: begin
    path ".command.err"     , emit: err
    path ".command.log"     , emit: log
    path ".command.out"     , emit: out
    path ".command.run"     , emit: run
    path ".command.sh"      , emit: sh
    path ".command.trace"   , emit: trace
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def query_reads = meta.single_end ? "${reads}" : "--first-pairs ${reads[0]} --second-pairs ${reads[1]}"
    """
    sylph \\
        profile \\
        $reference \\
        $query_reads \\
        -t ${task.cpus} \\
        ${task.ext.args} \\
        --output-file ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(echo \$(sylph --version 2>&1) | sed 's/^.*sylph //;s/ .*\$//')
    END_VERSIONS
    """
}
