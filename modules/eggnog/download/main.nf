process EGGNOG_DOWNLOAD {
    label 'process_low'
    label 'process_long'
    publishDir "${params.eggnog_db}", mode: params.publish_dir_mode, overwrite: true

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}":"${task.ext.docker}" }"

    output:
    path("eggnog/")     , emit: db, optional: true
    path("eggnog.tar.gz"), emit: db_tarball, optional: true
    path "*.{log,err}"   , emit: logs, optional: true
    path ".command.{begin,err,log,out,run,sh,trace}", emit: nf_logs
    path "versions.yml"  , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    echo "task.ext.args: ${task.ext.args}"
    
    # Create .command.begin
    date > .command.begin
    
    mkdir eggnog
    download_eggnog_data.py \\
        $args \\
        -y \\
        --data_dir eggnog/

    if [ "${task.ext.eggnog_save_as_tarball}" == "true" ]; then
        tar -czf eggnog.tar.gz eggnog/
        rm -rf eggnog/
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog-mapper: \$( echo \$(emapper.py --version 2>&1)| sed 's/.* emapper-//;s/ .*//')
    END_VERSIONS
    """
}
