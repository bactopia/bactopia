/**
 * Download and prepare database for Ariba analysis.
 *
 * This process executes ariba_getref to perform analysis
 *
 * @status stable
 * @keywords fastq, assembly, resistance, virulence
 * @tags complexity:moderate input-type:single output-type:single features:archive-output, compression, database-dependent, resource-download
 * @citation ariba_getref
 *
 * @note Requires external database to be available
 *
 * @input db_name
 * A database to setup up for Ariba
 *
 * @output db   An Ariba prepared database
 * @output logs Optional software execution logs containing warnings/errors
 */
nextflow.preview.types = true

process ARIBA_GETREF {
    tag "${db_name}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    db_name : String

    output:
    db   = file("ariba/ariba-${db_name}.tar.gz")
    logs = files("ariba/logs/*", optional: true)

    script:
    """
    # Download, format database, and tarball it
    ariba \\
        getref \\
        ${db_name} \\
        ${db_name}

    ariba \\
        prepareref \\
        -f ${db_name}.fa \\
        -m ${db_name}.tsv \\
        ${db_name}

    mv ${db_name}/ ariba-${db_name}/
    tar -zcvf ariba-${db_name}.tar.gz ariba-${db_name}/

    # Cleanup
    rm -rf ariba-${db_name}/

    # Move outputs to tool specific folder
    mkdir -p ariba/logs/${db_name}
    mv ariba-${db_name}.tar.gz ariba/
    cp .command.begin ariba/logs/${db_name}/nf.command.begin
    cp .command.err ariba/logs/${db_name}/nf.command.err
    cp .command.log ariba/logs/${db_name}/nf.command.log
    cp .command.out ariba/logs/${db_name}/nf.command.out
    cp .command.run ariba/logs/${db_name}/nf.command.run
    cp .command.sh ariba/logs/${db_name}/nf.command.sh
    cp .command.trace ariba/logs/${db_name}/nf.command.trace

    cat <<-END_VERSIONS > ariba/logs/${db_name}/versions.yml
    "${task.process}":
        ariba:  \$(echo \$(ariba version 2>&1) | sed 's/^.*ARIBA version: //;s/ .*\$//')
    END_VERSIONS
    """
}
