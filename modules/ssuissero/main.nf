process SSUISSERO {
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
    def VERSION = '1.0.1' // Version information not provided by tool on CLI
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    SsuisSero.sh \\
        -i $fasta_name \\
        -o ./ \\
        -s $prefix \\
        -x fasta \\
        -t $task.cpus

    # Cleanup
    rm -rf ${fasta_name} blast_res/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ssuissero: $VERSION
    END_VERSIONS
    """
}
