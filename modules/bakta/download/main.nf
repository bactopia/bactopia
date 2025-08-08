process BAKTA_DOWNLOAD {
    label 'process_low'
    label 'process_long'

    publishDir "${params.bakta_db}", mode: params.publish_dir_mode, overwrite: true
    
    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    output:
    path "bakta-${params.bakta_db_type}/*"     , emit: db, optional: true
    path "bakta-${params.bakta_db_type}.tar.gz", emit: db_tarball, optional: true
    path "*.{log,err}"                          , emit: logs, optional: true
    path ".command.begin"                       , emit: nf_begin
    path ".command.err"                         , emit: nf_err
    path ".command.log"                         , emit: nf_log
    path ".command.out"                         , emit: nf_out
    path ".command.run"                         , emit: nf_run
    path ".command.sh"                          , emit: nf_sh
    path ".command.trace"                       , emit: nf_trace
    path "versions.yml"                         , emit: versions

    script:
    """
    bakta_db \\
        download \\
        --type "${params.bakta_db_type}" \\
        --output bakta

    if [ "${params.bakta_save_as_tarball}" == "true" ]; then
        tar -czf bakta-${params.bakta_db_type}.tar.gz bakta/
        rm -rf bakta/
    else
        mv bakta/ bakta-${params.bakta_db_type}/
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$( echo \$(bakta --version 2>&1) | sed 's/^.*bakta //' )
    END_VERSIONS
    """
}
