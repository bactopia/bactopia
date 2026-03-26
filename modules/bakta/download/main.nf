/**
 * Download the Bakta annotation database.
 *
 * Fetches the pre-compiled database required by [Bakta](https://github.com/oschwengers/bakta)
 * for genome annotation. The database contains UniProt clusters, AMR genes, and other
 * reference data needed for comprehensive bacterial genome annotation.
 *
 * @status stable
 * @keywords bacteria, database, download, annotation, bakta, setup
 * @tags complexity:simple input-type:none output-type:multiple features:internet-access,archive-output,compression,resource-download,no-test
 * @citation bakta
 *
 * @note Internet & Storage Required
 * This process requires an active internet connection and significant disk space
 * to store the database files. The 'light' database is ~1.5GB, while the 'full'
 * database is ~30GB uncompressed.
 *
 * @output record(db, db_tarball, logs)
 * - `db`: The Bakta database directory containing annotation reference data
 * - `db_tarball`: A compressed tarball of the database (if requested via parameters)
 */
nextflow.preview.types = true

// bactopia-lint: ignore M012
process BAKTA_DOWNLOAD {
    label 'process_low'
    label 'process_long'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    output:
    record(
        db:         file("bakta-${task.ext.bakta_db_type}", optional: true),
        db_tarball: file("bakta-${task.ext.bakta_db_type}.tar.gz", optional: true),
        logs:       files("logs/*", optional: true)
    )

    script:
    """
    bakta_db \\
        download \\
        --type "${task.ext.bakta_db_type}" \\
        --output bakta

    if [ "${task.ext.save_as_tarball}" == "true" ]; then
        tar -czf bakta-${task.ext.bakta_db_type}.tar.gz bakta/
        rm -rf bakta/
    else
        mv bakta/ bakta-${task.ext.bakta_db_type}/
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
        bakta: \$( echo \$(bakta --version 2>&1) | sed 's/^.*bakta //' )
    END_VERSIONS
    """
}
