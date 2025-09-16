process BAKTA_DOWNLOAD {
    label 'process_low'
    label 'process_long'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    output:
    path "bakta-${task.ext.bakta_db_type}/*", emit: db, optional: true
    path "bakta-${task.ext.bakta_db_type}.tar.gz", emit: db_tarball, optional: true
    path "*.{log,err}", emit: logs, optional: true
    path ".command.begin", emit: begin, optional: true
    path ".command.err", emit: err
    path ".command.log", emit: log_file
    path ".command.out", emit: out
    path ".command.run", emit: run
    path ".command.sh", emit: sh
    path ".command.trace", emit: trace
    path "versions.yml", emit: versions

    script:
    """
    bakta_db \\
        download \\
        --type "${task.ext.bakta_db_type}" \\
        --output bakta

    if [ "${task.ext.save_as_tarball}" == "true" ]; then
        tar -czf bakta-${task.ext.bakta_db_type}.tar.gz bakta/
        rm -rf bakta/
    else
        mv bakta/ bakta-${task.ext.bakta_db_type}/
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$( echo \$(bakta --version 2>&1) | sed 's/^.*bakta //' )
    END_VERSIONS
    """
}
