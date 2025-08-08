VERSION = '2.2.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
process SRAHUMANSCRUBBER_INITDB {
    label 'process_single'
    storeDir params.datasets_cache
    publishDir params.datasets_cache
    
    conda "${task.ext.conda_env}"
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-human-scrubber:2.2.1--hdfd78af_0' :
        'quay.io/biocontainers/sra-human-scrubber:2.2.1--hdfd78af_0' }"

    output:
    path "*.human_filter.db", emit: db
    path "*.{log,err}"      , emit: logs, optional: true
    path ".command.begin"   , emit: begin
    path ".command.err"     , emit: err
    path ".command.log"     , emit: log
    path ".command.out"     , emit: out
    path ".command.run"     , emit: run
    path ".command.sh"      , emit: sh
    path ".command.trace"   , emit: trace
    path "versions.yml"     , emit: versions

    script:
    """
    echo "DOWNLOAD_START" > .command.begin
    DBVERSION=\$(curl "https://ftp.ncbi.nlm.nih.gov/sra/dbs/human_filter/current/version.txt")
    curl -f "https://ftp.ncbi.nlm.nih.gov/sra/dbs/human_filter/\${DBVERSION}.human_filter.db" -o "\${DBVERSION}.human_filter.db"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sra-human-scrubber: $VERSION
        sra-human-scrubber-db: \$DBVERSION
    END_VERSIONS
    """
}
