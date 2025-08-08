process SCCMEC {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}" :
        "${task.ext.docker}" }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    tuple val(meta), path("*.targets.blastn.tsv") , emit: targets
    tuple val(meta), path("*.targets.details.tsv"), emit: target_details
    tuple val(meta), path("*.regions.blastn.tsv") , emit: regions
    tuple val(meta), path("*.regions.details.tsv"), emit: regions_details
    path "versions.yml", emit: versions
    path ".command.begin", emit: begin
    path ".command.err", emit: err
    path ".command.log", emit: log
    path ".command.out", emit: out
    path ".command.run", emit: run
    path ".command.sh", emit: sh
    path ".command.trace", emit: trace

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    sccmec \\
        $args \\
        --prefix $prefix \\
        --input $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sccmec: \$(echo \$(sccmec --version 2>&1) | sed 's/.*sccmec_regions, version //;s/ .*\$//' )
        camlhmp: \$(echo \$(sccmec --version 2>&1) | sed 's/.*camlhmp, version //;s/ schema.*\$//' )
    END_VERSIONS
    """
}
