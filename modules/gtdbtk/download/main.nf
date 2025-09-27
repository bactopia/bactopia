process GTDBTK_DOWNLOAD {
    label 'process_low'
    label 'process_long'
    publishDir "${params.gtdb}", mode: params.publish_dir_mode, overwrite: true

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    output:
    path "gtdbtk/"      , emit: db, optional: true
    path "gtdbtk.tar.gz", emit: db_tarball, optional: true
    path "logs/*"       , emit: logs, optional: true

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
