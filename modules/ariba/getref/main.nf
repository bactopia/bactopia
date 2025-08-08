process ARIBA_GETREF {
    tag "$db_name"
    label 'process_low'
    storeDir params.datasets_cache
    publishDir params.datasets_cache

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    val(db_name)

    output:
    path("ariba-${db_name}.tar.gz") , emit: db
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
    # Download, format database, and tarball it
    ariba \\
        getref \\
        ${db_name} \\
        ${db_name}

    ariba \\
        prepareref \\
        -f ${db_name}.fa \\
        -m ${db_name}.tsv \\
        ${db_name}

    mv ${db_name}/ ariba-${db_name}/
    tar -zcvf ariba-${db_name}.tar.gz ariba-${db_name}/

    # Cleanup
    rm -rf ariba-${db_name}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ariba:  \$(echo \$(ariba version 2>&1) | sed 's/^.*ARIBA version: //;s/ .*\$//')
    END_VERSIONS
    """
}
