/**
 * Download the eggNOG database for functional annotation.
 *
 * Fetches the pre-computed orthology data and Diamond database required by
 * [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper). This includes the massive
 * protein database and taxonomic information needed for accurate ortholog assignment.
 *
 * @status stable
 * @keywords eggnog, database, download, annotation, functional, orthology
 * @tags complexity:simple input-type:none output-type:multiple features:internet-access,resource-download,no-test
 * @citation eggnog_mapper
 *
 * @note Internet & Storage Required
 * This process requires an active internet connection and significant disk space (often >50GB)
 * to store the uncompressed database files.
 *
 * @output record(db, db_tarball, logs)
 * - `db`: The eggNOG database directory (Diamond database and taxonomy info)
 * - `db_tarball`: A compressed tarball of the database (if requested via parameters)
 */
nextflow.preview.types = true

// bactopia-lint: ignore M012,M017,M018,M022,M023,M024,M025,M026,M028
process EGGNOG_DOWNLOAD {
    label 'process_low'
    label 'process_long'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    output:
    record(
        db: file("eggnog", optional: true),
        db_tarball: file("eggnog.tar.gz", optional: true),
        logs: files("logs/*", optional: true)
    )

    script:
    def args = task.ext.args ?: ''
    """
    mkdir eggnog
    download_eggnog_data.py \\
        ${args} \\
        -y \\
        --data_dir eggnog/

    if [ "${task.ext.eggnog_save_as_tarball}" == "true" ]; then
        tar -czf eggnog.tar.gz eggnog/
        rm -rf eggnog/
    fi

    # Capture logs
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
        eggnog-mapper: \$( echo \$(emapper.py --version 2>&1)| sed 's/.* emapper-//;s/ .*//')
    END_VERSIONS
    """
}
