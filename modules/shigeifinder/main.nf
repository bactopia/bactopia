process SHIGEIFINDER {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}" :
        "${task.ext.docker}" }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions
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
    def VERSION = '1.3.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    shigeifinder \\
        $args \\
        --output ${prefix}.tsv \\
        -t $task.cpus \\
        -i $fasta_name

    # Cleanup
    rm -rf ${fasta_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shigeifinder: $VERSION
    END_VERSIONS
    """
}
