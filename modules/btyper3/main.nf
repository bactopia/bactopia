process BTYPER3 {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("results/*_final_results.txt"), emit: tsv
    tuple val(meta), path("results/*")                  , emit: results
    path "*.{log,err}" , emit: logs, optional: true
    tuple val(meta), path(".command.begin")                              , emit: nf_begin
    tuple val(meta), path(".command.err")                                , emit: nf_err
    tuple val(meta), path(".command.log")                                , emit: nf_log
    tuple val(meta), path(".command.out")                                , emit: nf_out
    tuple val(meta), path(".command.run")                                , emit: nf_run
    tuple val(meta), path(".command.sh")                                 , emit: nf_sh
    tuple val(meta), path(".command.trace")                              , emit: nf_trace
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    meta.output_dir = "${meta.id}/main/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/main/${task.ext.process_name}/${task.ext.subdir}/logs"
    meta.process_name = task.ext.process_name
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    # Btyper3 does not accept compressed files
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    btyper3 \\
        $task.ext.args \\
        --output ./ \\
        --input ${fasta_name}

    mv btyper3_final_results/ results/

    # Cleanup
    rm -rf ${fasta_name} ${fasta_name}.njs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        btyper3: \$(echo \$(btyper3 --version 2>&1) | sed 's/^.*btyper3 //;' ))
    END_VERSIONS
    """
}
