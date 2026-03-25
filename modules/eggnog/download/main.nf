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
 * @output record(db, logs)
 * - `db`: The eggNOG database files (Diamond database and taxonomy info)
 */
nextflow.preview.types = true

// bactopia-lint: ignore M012
process EGGNOG_DOWNLOAD {
    label 'process_low'
    label 'process_long'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    output:
    record(
        db:   files("eggnog/eggnog*"),
        logs: files("eggnog/logs/*")
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
        mkdir eggnog
        mv eggnog.tar.gz eggnog/
    fi

    # Capture logs
    mkdir eggnog/logs
    cp .command.begin eggnog/logs/nf.command.begin
    cp .command.err eggnog/logs/nf.command.err
    cp .command.log eggnog/logs/nf.command.log
    cp .command.out eggnog/logs/nf.command.out
    cp .command.run eggnog/logs/nf.command.run
    cp .command.sh eggnog/logs/nf.command.sh
    cp .command.trace eggnog/logs/nf.command.trace

    cat <<-END_VERSIONS > eggnog/logs/versions.yml
    "${task.process}":
        eggnog-mapper: \$( echo \$(emapper.py --version 2>&1)| sed 's/.* emapper-//;s/ .*//')
    END_VERSIONS
    """
}
