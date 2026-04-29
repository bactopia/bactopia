/**
 * Download and configure the GTDB-Tk reference database.
 *
 * Uses the official `download-db.sh` script to fetch the latest Genome Taxonomy Database (GTDB)
 * files required by [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk). It automatically uncompresses
 * the data and verifies the installation using `gtdbtk check_install`.
 *
 * @status stable
 * @keywords gtdb, taxonomy, database, download, setup, bacteria, archaea
 * @tags complexity:simple input-type:none output-type:multiple features:internet-access,resource-download,conditional-logic,no-test
 * @citation gtdb_tk
 *
 * @note Internet & Storage Required
 * This process requires an active internet connection and significant disk space (~60GB+ uncompressed)
 * to store the database files.
 *
 * @output record(db, db_tarball, logs)
 * - `db`: The directory containing the uncompressed GTDB-Tk database files
 * - `db_tarball`: A compressed tarball of the database (if requested via parameters)
 */
nextflow.enable.types = true

// bactopia-lint: ignore M012,M017,M018,M022,M023,M024,M025,M026,M028,M033
process GTDBTK_DOWNLOAD {
    label 'process_low'
    label 'process_long'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    output:
    record(
        db: file("gtdbtk", optional: true),
        db_tarball: file("gtdbtk.tar.gz", optional: true),
        logs: files("logs/*", optional: true)
    )

    script:
    """
    # Create .command.begin
    date > .command.begin

    export GTDBTK_DATA_PATH="./gtdbtk"
    mkdir ./gtdbtk
    download-db.sh ./gtdbtk
    gtdbtk check_install && touch gtdb-setup.txt

    if [ "${task.ext.gtdb_save_as_tarball}" == "true" ]; then
        tar -czf gtdbtk.tar.gz gtdbtk/
        rm -rf gtdbtk/
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
        gtdbtk: \$(echo \$(gtdbtk --version -v 2>&1) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """
}
