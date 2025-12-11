/**
 * Download and setup GTDB-Tk database.
 *
 * This process executes gtdbtk_download to perform analysis
 *
 * @status stable
 * @keywords gtdb, database, download, setup
 * @tags complexity:moderate input-type:single output-type:multiple features:archive-output, compression, conditional-logic, resource-download
 * @citation gtdbtk_download
 * @output db         GTDB database directory
 * @output db_tarball GTDB database tarball (optional)
 * @output logs       Optional tool execution logs
 */
nextflow.preview.types = true

process GTDBTK_DOWNLOAD {
    label 'process_low'
    label 'process_long'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    output:
    db         = files("gtdbtk/", optional: true)
    db_tarball = files("gtdbtk.tar.gz", optional: true)
    logs       = files("logs/*", optional: true)

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
