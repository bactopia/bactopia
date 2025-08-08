process HPSUISSERO {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}":"${task.ext.docker}" }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "*.{log,err}"            , emit: logs, optional: true
    path ".command.{begin,err,log,out,run,sh,trace}", emit: nf_logs
    path "versions.yml"           , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    echo "task.ext.args: ${task.ext.args}"
    
    # Create .command.begin
    date > .command.begin
    
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    HpsuisSero.sh \\
        -i $fasta_name \\
        -o ./ \\
        -s $prefix \\
        -x fasta \\
        -t $task.cpus

    # Cleanup
    rm -rf ${fasta_name} blast_res/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hpsuissero: ${task.ext.version}
    END_VERSIONS
    """
}
