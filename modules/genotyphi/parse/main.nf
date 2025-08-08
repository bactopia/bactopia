process GENOTYPHI_PARSE {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}":"${task.ext.docker}" }"

    input:
    tuple val(meta), path(json)

    output:
    tuple val(meta), path("*.tsv")     , emit: tsv
    path "*.{log,err}", optional: true , emit: logs
    path ".command.{begin,err,log,out,run,sh,trace}", emit: nf_logs
    path "versions.yml"                , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "task.ext.args: ${task.ext.args}"
    
    # Create .command.begin
    date > .command.begin
    
    parse_typhi_mykrobe.py \\
        --jsons $json \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genotyphi: \$(echo \$(genotyphi --version 2>&1) | sed 's/^.*GenoTyphi v//;' )
    END_VERSIONS
    """
}
