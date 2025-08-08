process KLEBORATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}":"${task.ext.docker}" }"

    input:
    tuple val(meta), path(fastas)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "*.{log,err}", emit: logs, optional: true
    path ".command.{begin,err,log,out,run,sh,trace}", emit: nf_logs
    path "versions.yml",emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "task.ext.args: ${task.ext.args}"
    
    # Create .command.begin
    date > .command.begin
    
    mkdir results/
    kleborate \\
        $args \\
        --outdir results/ \\
        --assemblies $fastas

    # Rename output file to include the prefix name
    find results/ -name "*output.txt" -print0 | while read -d \$'\0' file; do mv "\$file" "${prefix}.txt"; done

    # Negative results will not have an output file
    if [ ! -f "${prefix}.txt" ]; then
        touch "${prefix}.txt"
    fi

    # cleanup
    rm -rf results/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kleborate: \$( echo \$(kleborate --version | sed 's/Kleborate v//;'))
    END_VERSIONS
    """
}
