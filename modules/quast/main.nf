process QUAST {
    tag "$meta.id"
    label 'process_medium'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(fasta), path(meta_file)

    output:
    tuple val(meta), path("results/${meta.id}.tsv"), emit: tsv
    tuple val(meta), path("results/*")             , emit: results
    tuple val(meta), path("*.{log,err}")           , emit: logs, optional: true
    tuple val(meta), path(".command.begin")        , emit: nf_begin
    tuple val(meta), path(".command.err")          , emit: nf_err
    tuple val(meta), path(".command.log")          , emit: nf_log
    tuple val(meta), path(".command.out")          , emit: nf_out
    tuple val(meta), path(".command.run")          , emit: nf_run
    tuple val(meta), path(".command.sh")           , emit: nf_sh
    tuple val(meta), path(".command.trace")        , emit: nf_trace
    tuple val(meta), path("versions.yml")          , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs"
    meta.process_name = task.ext.process_name
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
