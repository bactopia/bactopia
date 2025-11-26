nextflow.preview.types = true

process DEFENSEFINDER_UPDATE {
    tag "update"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    output:
    db   = file("defense-finder/defense-finder-models-${task.ext.df_models_version}.tar")
    logs = files("defense-finder/logs/*", optional: true)

    script:
    """
    mkdir models
    wget \\
        -O models/defense-finder-models-v${task.ext.df_models_version}.tar.gz \\
        https://github.com/mdmparis/defense-finder-models/archive/refs/tags/${task.ext.df_models_version}.tar.gz

    wget \\
        -O models/CasFinder-${task.ext.casfinder_version}.tar.gz \\
        https://github.com/macsy-models/CasFinder/archive/refs/tags/${task.ext.casfinder_version}.tar.gz

    tar -cvf defense-finder-models-${task.ext.df_models_version}.tar models/

    # Move outputs to tool specific folder
    mkdir -p defense-finder/logs
    mv defense-finder-models-${task.ext.df_models_version}.tar defense-finder/
    cp .command.begin defense-finder/logs/nf.command.begin
    cp .command.err defense-finder/logs/nf.command.err
    cp .command.log defense-finder/logs/nf.command.log
    cp .command.out defense-finder/logs/nf.command.out
    cp .command.run defense-finder/logs/nf.command.run
    cp .command.sh defense-finder/logs/nf.command.sh
    cp .command.trace defense-finder/logs/nf.command.trace

    cat <<-END_VERSIONS > defense-finder/logs/versions.yml
    "${task.process}":
        defense-finder: ${task.ext.df_version}
        defense-finder-models: ${task.ext.df_models_version}
        casfinder-models: ${task.ext.casfinder_version}
    END_VERSIONS
    """
}
