/**
 * Download eggNOG database files for annotation.
 *
 * This process executes eggnog_download to perform analysis
 *
 * @status stable
 * @keywords eggnog, database, download, annotation
 * @tags complexity:moderate input-type:single output-type:single features:archive-output, compression, conditional-logic, resource-download
 * @citation eggnog_download
 * @output db   Directory containing downloaded eggNOG database files
 * @output logs Optional tool execution logs
 */
nextflow.preview.types = true

process EGGNOG_DOWNLOAD {
    label 'process_low'
    label 'process_long'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    output:
    db   = files("eggnog/eggnog*")
    logs = files("eggnog/logs/*")

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
