/**
 * Download the nohuman database for human read removal.
 *
 * Fetches the Kraken2-based database used by [nohuman](https://github.com/mbhall88/nohuman)
 * to classify and remove human reads from sequencing datasets. The database is built from
 * Human Pangenome Reference Consortium (HPRC) genomes.
 *
 * @status stable
 * @keywords bacteria, database, download, human, decontamination, kraken2, nohuman
 * @tags complexity:simple input-type:none output-type:multiple features:internet-access,archive-output,compression,resource-download,no-test
 * @citation kraken2
 *
 * @note Internet & Storage Required
 * This process requires an active internet connection and sufficient disk space
 * to store the Kraken2 database files.
 *
 * @output record(db, db_tarball, logs)
 * - `db`: The nohuman Kraken2 database directory
 * - `db_tarball`: A compressed tarball of the database (if requested via parameters)
 */
nextflow.enable.types = true

// bactopia-lint: ignore M012,M017,M018,M022,M023,M024,M025,M026,M028,M033
process NOHUMAN_DOWNLOAD {
    label 'process_low'
    label 'process_long'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    output:
    record(
        db: file("nohuman-db", optional: true),
        db_tarball: file("nohuman-db.tar.gz", optional: true),
        logs: files("logs/*", optional: true)
    )

    script:
    """
    nohuman \\
        --download \\
        --db nohuman-db \\
        ${task.ext.args}

    if [ "${task.ext.nohuman_save_as_tarball}" == "true" ]; then
        tar -czf nohuman-db.tar.gz nohuman-db/
        rm -rf nohuman-db/
    fi

    # Move outputs to tool specific folder
    mkdir -p logs
    cp .command.begin logs/nf.command.begin
    cp .command.err logs/nf.command.err
    cp .command.log logs/nf.command.log
    cp .command.out logs/nf.command.out
    cp .command.run logs/nf.command.run
    cp .command.sh logs/nf.command.sh
    cp .command.trace logs/nf.command.trace

    cat <<-END_VERSIONS > logs/versions.yml
    "${task.process}":
        nohuman: \$( nohuman --version 2>&1 | sed 's/nohuman //' )
    END_VERSIONS
    """
}
