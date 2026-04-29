/**
 * Download the MIDAS reference database.
 *
 * Fetches the pre-compiled database required by [MIDAS](https://github.com/snayfach/MIDAS)
 * for metagenomic species profiling. The database contains reference genomes from the
 * UHGG collection used for species identification and abundance estimation.
 *
 * @status stable
 * @keywords midas, download, database, metagenomics, species
 * @tags complexity:simple input-type:none output-type:multiple features:internet-access,resource-download,no-test
 * @citation midas
 *
 * @note Internet & Storage Required
 * This process requires an active internet connection and significant disk space
 * to store the database files (~3.3GB compressed, ~7GB uncompressed).
 *
 * @output record(db, db_tarball, logs)
 * - `db`: The MIDAS database directory containing reference genome data
 * - `db_tarball`: A compressed tarball of the database (if requested via parameters)
 */
nextflow.enable.types = true

// bactopia-lint: ignore M012,M017,M018,M022,M023,M024,M025,M026,M028,M033
process MIDAS_DOWNLOAD {
    label 'process_low'
    label 'process_long'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    output:
    record(
        db: file("midas_db_v${db_version}", optional: true),
        db_tarball: file("midas_db_v${db_version}.tar.gz", optional: true),
        logs: files("logs/*", optional: true)
    )

    script:
    db_version = "1.2"
    """
    wget ${task.ext.args} -O midas_db_v${db_version}.tar.gz \
        https://midasdb.pollard.gladstone.org/uhgg/midas_db_v${db_version}.tar.gz

    if [ "${task.ext.save_as_tarball}" == "true" ]; then
        echo "Keeping database as tarball"
    else
        tar -xzf midas_db_v${db_version}.tar.gz
        rm -f midas_db_v${db_version}.tar.gz
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
        wget: \$(echo \$( wget --version 2>&1) | sed "s/GNU Wget //;s/ .*//;" )
        midas_db: ${db_version}
    END_VERSIONS
    """
}
