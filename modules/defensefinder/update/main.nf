/**
 * Download and package the DefenseFinder and CasFinder model databases.
 *
 * Fetches the latest HMM profiles from the [DefenseFinder](https://github.com/mdmparis/defense-finder-models)
 * and [CasFinder](https://github.com/macsy-models/CasFinder) repositories, then packages them
 * into a single tarball for use by the DefenseFinder module.
 *
 * @status stable
 * @keywords bacteria, defense, database, download, hmm, casfinder, crispr
 * @tags complexity:simple input-type:none output-type:single features:internet-access,resource-download,archive-output
 * @citation defensefinder
 *
 * @note Internet Required
 * This process requires an active internet connection to fetch the models from GitHub.
 *
 * @output record(db, logs)
 * - `db`: A compressed tarball containing both DefenseFinder and CasFinder models
 */
nextflow.preview.types = true

// bactopia-lint: ignore M012,M017,M018,M023,M024,M025,M026,M028
process DEFENSEFINDER_UPDATE {
    tag "update"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    output:
    record(
        db: file("defense-finder/defense-finder-models-${task.ext.defensefinder_models_version}.tar"),
        logs: files("defense-finder/logs/*", optional: true)
    )

    script:
    """
    mkdir models
    wget \\
        -O models/defense-finder-models-v${task.ext.defensefinder_models_version}.tar.gz \\
        https://github.com/mdmparis/defense-finder-models/archive/refs/tags/${task.ext.defensefinder_models_version}.tar.gz

    wget \\
        -O models/CasFinder-${task.ext.casfinder_version}.tar.gz \\
        https://github.com/macsy-models/CasFinder/archive/refs/tags/${task.ext.casfinder_version}.tar.gz

    tar -cvf defense-finder-models-${task.ext.defensefinder_models_version}.tar models/

    # Move outputs to tool specific folder
    mkdir -p defense-finder/logs
    mv defense-finder-models-${task.ext.defensefinder_models_version}.tar defense-finder/
    cp .command.begin defense-finder/logs/nf.command.begin
    cp .command.err defense-finder/logs/nf.command.err
    cp .command.log defense-finder/logs/nf.command.log
    cp .command.out defense-finder/logs/nf.command.out
    cp .command.run defense-finder/logs/nf.command.run
    cp .command.sh defense-finder/logs/nf.command.sh
    cp .command.trace defense-finder/logs/nf.command.trace

    # Cleanup
    rm -rf models/

    cat <<-END_VERSIONS > defense-finder/logs/versions.yml
    "${task.process}":
        defense-finder: ${task.ext.defensefinder_version}
        defense-finder-models: ${task.ext.defensefinder_models_version}
        casfinder-models: ${task.ext.casfinder_version}
    END_VERSIONS
    """
}
