process DEFENSEFINDER_UPDATE {
    tag "update"
    label 'process_low'
    storeDir params.datasets_cache
    publishDir params.datasets_cache

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}":"${task.ext.docker}" }"

    output:
    path "defense-finder-models-${task.ext.df_models_version}.tar", emit: db

    script:
    """
    echo "task.ext.args: ${task.ext.args}"
    
    # Create .command.begin
    date > .command.begin
    
    mkdir models
    wget \\
        -O models/defense-finder-models-v${task.ext.df_models_version}.tar.gz \\
        https://github.com/mdmparis/defense-finder-models/archive/refs/tags/${task.ext.df_models_version}.tar.gz

    wget \\
        -O models/CasFinder-${task.ext.casfinder_version}.tar.gz \\
        https://github.com/macsy-models/CasFinder/archive/refs/tags/${task.ext.casfinder_version}.tar.gz

    tar -cvf defense-finder-models-${task.ext.df_models_version}.tar models/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        defense-finder: ${task.ext.df_version}
        defense-finder-models: ${task.ext.df_models_version}
        casfinder-models: ${task.ext.casfinder_version}
    END_VERSIONS
    """
}
