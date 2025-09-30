process AMRFINDERPLUS_UPDATE {
    tag "update"
    label 'process_low'

    conda "${task.ext.conda_env}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        task.ext.container_image : 
        task.ext.container }"

    output:
    path "amrfinderplus.tar.gz", emit: db
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
    mkdir amrfinderplus-temp
    amrfinder_update -d amrfinderplus-temp
    mv amrfinderplus-temp/\$(readlink amrfinderplus-temp/latest) amrfinderplus/
    tar czvf amrfinderplus.tar.gz amrfinderplus/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
        amrfinderplus-database: \$(echo \$(echo \$(amrfinder --database amrfinderplus --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev))
    END_VERSIONS
    """
}
