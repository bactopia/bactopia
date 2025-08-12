process CHECKM2_DOWNLOAD {
    label 'process_low'
    label 'process_long'

    publishDir "${params.checkm2_db}", mode: params.publish_dir_mode, overwrite: true

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    output:
    path "checkm2_db_v${db_version}.dmnd"  , emit: db
    path "contents.json"                   , emit: json
    tuple val(meta), path("*.{log,err}")   , emit: logs, optional: true
    tuple val(meta), path(".command.begin"), emit: nf_begin
    tuple val(meta), path(".command.err")  , emit: nf_err
    tuple val(meta), path(".command.log")  , emit: nf_log
    tuple val(meta), path(".command.out")  , emit: nf_out
    tuple val(meta), path(".command.run")  , emit: nf_run
    tuple val(meta), path(".command.sh")   , emit: nf_sh
    tuple val(meta), path(".command.trace"), emit: nf_trace
    tuple val(meta), path("versions.yml")  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs"
    meta.process_name = task.ext.process_name
    zenodo_id  = 5571251  // Default to latest version 
    api_data   = (new groovy.json.JsonSlurper()).parseText(file("https://zenodo.org/api/records/${zenodo_id}").text)
    db_version = api_data.metadata.version
    checksum   = api_data.files[0].checksum.replaceFirst(/^md5:/, "md5=")
    meta       = [id: 'checkm2_db', version: db_version]
    """
    # Automatic download is broken when using singularity/apptainer (https://github.com/chklovski/CheckM2/issues/73)
    # So it's necessary to download the database manually
    aria2c \\
        $task.ext.args \\
        --checksum ${checksum} \\
        https://zenodo.org/records/${zenodo_id}/files/checkm2_database.tar.gz

    echo ${db_version}

    tar -xzf checkm2_database.tar.gz
    db_path=\$(find -name *.dmnd)
    mv \$db_path checkm2_db_v${db_version}.dmnd
    rmdir CheckM2_database

    mv CONTENTS.json contents.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aria2: \$(echo \$(aria2c --version 2>&1) | grep 'aria2 version' | cut -f3 -d ' ')
        checkm2_db: ${db_version}
    END_VERSIONS
    """

    stub:
    """
    touch checkm_db.dmnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
    END_VERSIONS
    """
}
