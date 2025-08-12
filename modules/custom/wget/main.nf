process CUSTOM_WGET {
    label 'process_low'
    storeDir params.datasets_cache
    publishDir params.datasets_cache

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    output:
    path "${task.ext.args2}", emit: download
    path "*.{log,err}"      , emit: logs, optional: true
    path ".command.begin"   , emit: nf_nf_begin
    path ".command.err"     , emit: nf_nf_err
    path ".command.log"     , emit: nf_nf_log
    path ".command.out"     , emit: nf_nf_out
    path ".command.run"     , emit: nf_nf_run
    path ".command.sh"      , emit: nf_nf_sh
    path ".command.trace"   , emit: nf_nf_trace
    path "versions.yml"     , emit: versions

    script:
    """
    wget $task.ext.args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    wget: \$(echo \$( wget --version 2>&1) | sed "s/GNU Wget //;s/ .*//;" )
    END_VERSIONS
    """
}
