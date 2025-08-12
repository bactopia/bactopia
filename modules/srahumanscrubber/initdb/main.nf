process SRAHUMANSCRUBBER_INITDB {
    label 'process_single'
    storeDir params.datasets_cache
    publishDir params.datasets_cache

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    output:
    path "*.human_filter.db", emit: db
    path "*.{log,err}"   , emit: logs, optional: true
    path ".command.begin", emit: begin, optional: true
    path ".command.err"  , emit: err
    path ".command.log"  , emit: log_file
    path ".command.out"  , emit: out
    path ".command.run"  , emit: run
    path ".command.sh"   , emit: sh
    path ".command.trace", emit: trace
    path "versions.yml"  , emit: versions

    script:
    """
    echo "DOWNLOAD_START" > .command.begin
    DBVERSION=\$(curl "https://ftp.ncbi.nlm.nih.gov/sra/dbs/human_filter/current/version.txt")
    curl -f "https://ftp.ncbi.nlm.nih.gov/sra/dbs/human_filter/\${DBVERSION}.human_filter.db" -o "\${DBVERSION}.human_filter.db"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    sra-human-scrubber: 2.2.1
    sra-human-scrubber-db: \$DBVERSION
    END_VERSIONS
    """
}
