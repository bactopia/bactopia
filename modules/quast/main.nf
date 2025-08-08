process QUAST {
    tag "$meta.id"
    label 'process_medium'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}" :
        "${task.ext.docker}" }"

    input:
    tuple val(meta), path(fasta), path(meta_file)

    output:
    tuple val(meta), path("results/${meta.id}.tsv"), emit: tsv
    tuple val(meta), path("results/*")             , emit: results
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
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    est_ref_size=""
    # Use rev to get the last column easily, then re-reverse it
    ref_size=\$(tail -n 1 $meta_file | rev | cut -f 1 | rev)
    if [ "\${ref_size}" != "0" ]; then
        est_ref_size="--est-ref-size \${ref_size}"
    fi

    quast ${fasta_name} \${est_ref_size} \\
        -o results \\
        --threads ${task.cpus} \\
        $args \\
        --glimmer

    mv results/quast.log ./
    mv results/transposed_report.tsv results/${prefix}.tsv

    # Cleanup
    rm -rf ${fasta_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
