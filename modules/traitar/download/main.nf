/**
 * Download the Pfam database required by Traitar.
 *
 * Fetches the pre-compiled Pfam database required by [Traitar](https://github.com/nick-youngblut/traitar3/)
 * for microbial phenotype prediction. The database contains Pfam HMM models used for
 * protein family annotation during trait prediction.
 *
 * @status stable
 * @keywords phenotype, traits, pfam, database, download
 * @tags complexity:simple input-type:none output-type:single features:internet-access,resource-download,no-test
 * @citation traitar
 *
 * @note Internet & Storage Required
 * This process requires an active internet connection and significant disk space
 * to store the Pfam database files (~1.2GB).
 *
 * @output record(db, logs)
 * - `db`: The Pfam-A HMM file for Traitar
 */
nextflow.enable.types = true

// bactopia-lint: ignore M012,M017,M018,M023,M024,M025,M026,M028,M033
process TRAITAR_DOWNLOAD {
    label 'process_low'
    label 'process_long'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    output:
    record(
        db: file("Pfam-A.hmm"),
        logs: files("logs/*", optional: true)
    )

    script:
    """
    traitar pfam pfam_data
    mv pfam_data/Pfam-A.hmm Pfam-A.hmm

    # Move outputs to tool specific folder
    mkdir -p logs
    cp .command.begin logs/nf.command.begin
    cp .command.err logs/nf.command.err
    cp .command.log logs/nf.command.log
    cp .command.out logs/nf.command.out
    cp .command.run logs/nf.command.run
    cp .command.sh logs/nf.command.sh
    cp .command.trace logs/nf.command.trace

    # Cleanup

    cat <<-END_VERSIONS > logs/versions.yml
    "${task.process}":
        traitar: \$( traitar --version 2>&1 | tail -1 )
    END_VERSIONS
    """
}
