process IQTREE {
    tag "$prefix"
    label 'process_medium'
    label 'process_long'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}":"${task.ext.docker}" }"

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("${prefix}*")                         , emit: results
    tuple val(meta), path("${prefix}.treefile")                 , emit: phylogeny
    tuple val(meta), path(alignment), path("${prefix}.treefile"), emit: aln_tree
    path "*.{log,err}"                                          , emit: logs, optional: true
    path ".command.{begin,err,log,out,run,sh,trace}"            , emit: nf_logs
    path "versions.yml"                                         , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "task.ext.args: ${task.ext.args}"
    
    # Create .command.begin
    date > .command.begin
    
    iqtree \\
        $args \\
        -s $alignment \\
        -nt $task.cpus \\
        -ntmax $task.cpus \\
        -pre $prefix

    # Only gzip files if they exist
    if [[ -f "${prefix}.alninfo" ]]; then
        gzip ${prefix}.alninfo
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
    END_VERSIONS
    """
}
