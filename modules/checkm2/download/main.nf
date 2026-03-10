/**
 * Download the pre-trained CheckM2 database.
 *
 * Fetches the required Diamond database used by [CheckM2](https://github.com/chklovski/CheckM2)
 * for genome quality prediction. This database contains the machine learning model training data.
 *
 * @status stable
 * @keywords checkm2, download, database, diamond, machine learning
 * @tags complexity:simple input-type:none output-type:multiple features:internet-access,resource-download
 * @citation checkm2
 *
 * @note Internet Required
 * This process requires an active internet connection to fetch the database.
 *
 * @output record(db, json, logs)
 * - `db`: The CheckM2 Diamond database file (*.dmnd)
 * - `json`: Metadata file describing the database contents
 */
nextflow.preview.types = true

process CHECKM2_DOWNLOAD {
    label 'process_low'
    label 'process_long'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    output:
    db   = file("checkm2_db_v${db_version}.dmnd")
    json = file("checkm2_db_v${db_version}.json")
    logs = files("logs/*", optional: true)

    script:
    // Check for latest versions at https://doi.org/10.5281/zenodo.4626518
    zenodo_id = 14897628
    db_version = 3
    checksum = "07c10655620843b517d0df0c160d911f"
    """
    # Automatic download is broken when using singularity/apptainer (https://github.com/chklovski/CheckM2/issues/73)
    # So it's necessary to download the database manually
    aria2c \\
        ${task.ext.args} \\
        --checksum=md5=${checksum} \\
        https://zenodo.org/records/${zenodo_id}/files/checkm2_database.tar.gz

    echo ${db_version}
    tar -xzf checkm2_database.tar.gz
    db_path=\$(find -name *.dmnd)
    mv \$db_path checkm2_db_v${db_version}.dmnd

    # Cleanup
    rm -rf CheckM2_database/ checkm2_database.tar.gz
    mv CONTENTS.json checkm2_db_v${db_version}.json

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
        aria2: \$(echo \$(aria2c --version 2>&1) | grep 'aria2 version' | cut -f3 -d ' ')
        checkm2_db: ${db_version}
    END_VERSIONS
    """

    stub:
    """
    touch checkm_db.dmnd
    touch contents.json
    mkdir -p logs
    cat <<-END_VERSIONS > logs/versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
    END_VERSIONS
    """
}
