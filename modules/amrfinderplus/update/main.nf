process AMRFINDERPLUS_UPDATE {
    tag "update"
    label 'process_low'
    storeDir params.datasets_cache
    publishDir params.datasets_cache

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    output:
    path "amrfinderplus.tar.gz"     , emit: db
    path "*.{log,err}"              , emit: logs, optional: true
    path ".command.begin"           , emit: nf_begin
    path ".command.err"             , emit: nf_err
    path ".command.log"             , emit: nf_log
    path ".command.out"             , emit: nf_out
    path ".command.run"             , emit: nf_run
    path ".command.sh"              , emit: nf_sh
    path ".command.trace"           , emit: nf_trace
    path "versions.yml"             , emit: versions

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
