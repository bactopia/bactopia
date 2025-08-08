process GTDBTK_SETUPDB {
    label 'process_low'
    label 'process_long'
    publishDir "${params.gtdb}", mode: params.publish_dir_mode, overwrite: true

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}":"${task.ext.docker}" }"

    output:
    path("gtdbtk/")     , emit: db, optional: true
    path("gtdbtk.tar.gz"), emit: db_tarball, optional: true
    path "*.{log,err}"   , emit: logs, optional: true
    path ".command.{begin,err,log,out,run,sh,trace}", emit: nf_logs
    path "versions.yml"  , emit: versions

    script:
    """
    echo "task.ext.args: ${task.ext.args}"
    
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo \$(gtdbtk --version -v 2>&1) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """
}
