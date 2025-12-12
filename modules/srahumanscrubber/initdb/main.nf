/**
 * Initialize human read removal database for SRA Human Scrubber.
 *
 * Uses [SRA Human Scrubber](https://github.com/ncbi/sra-human-scrubber) to download and
 * initialize the necessary k-mer database required for scrubbing human reads from
 * sequencing data.
 *
 * @status stable
 * @keywords human, database, scrubber, ncbi, download
 * @tags complexity:simple input-type:none output-type:single features:internet-access,resource-download
 * @citation srahumanscrubber_initdb
 *
 * @output db   The initialized SRA Human Scrubber database files
 * @output logs Optional software execution logs containing warnings/errors
 */
nextflow.preview.types = true

process SRAHUMANSCRUBBER_INITDB {
    label 'process_single'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    output:
    db   = files("${prefix}/*human_filter.db*")
    logs = files("${prefix}/logs/*", optional: true)

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
