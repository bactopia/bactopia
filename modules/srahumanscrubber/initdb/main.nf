process SRAHUMANSCRUBBER_INITDB {
    label 'process_single'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    output:
    path "${prefix}/*human_filter.db*", emit: db
    path "${prefix}/logs/*"           , emit: logs, optional: true

    script:
    prefix = task.ext.process_name
    """
    mkdir -p ${prefix}/logs
    DBVERSION=\$(curl "https://ftp.ncbi.nlm.nih.gov/sra/dbs/human_filter/current/version.txt")
    curl -f "https://ftp.ncbi.nlm.nih.gov/sra/dbs/human_filter/\${DBVERSION}" -o "${prefix}/\${DBVERSION}"

    # Move outputs to tool specific folder
    cp .command.begin ${prefix}/logs/nf.command.begin
    cp .command.err ${prefix}/logs/nf.command.err
    cp .command.log ${prefix}/logs/nf.command.log
    cp .command.out ${prefix}/logs/nf.command.out
    cp .command.run ${prefix}/logs/nf.command.run
    cp .command.sh ${prefix}/logs/nf.command.sh
    cp .command.trace ${prefix}/logs/nf.command.trace

    cat <<-END_VERSIONS > ${prefix}/logs/versions.yml
    "${task.process}":
        sra-human-scrubber: 2.2.1
        sra-human-scrubber-db: \$DBVERSION
    END_VERSIONS
    """
}
