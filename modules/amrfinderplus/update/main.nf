nextflow.preview.types = true

process AMRFINDERPLUS_UPDATE {
    tag "amrfinderplus-update"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    output:
    db   = file("updater/amrfinderplus.tar.gz")
    logs = files("updater/logs/*", optional: true)

    script:
    """
    mkdir -p updater/logs
    mkdir amrfinderplus-temp
    amrfinder_update -d amrfinderplus-temp
    mv amrfinderplus-temp/\$(readlink amrfinderplus-temp/latest) amrfinderplus/
    tar czvf amrfinderplus.tar.gz amrfinderplus/
    mv amrfinderplus.tar.gz updater/ 

    # Move outputs to tool specific folder
    cp .command.begin updater/logs/nf.command.begin
    cp .command.err updater/logs/nf.command.err
    cp .command.log updater/logs/nf.command.log
    cp .command.out updater/logs/nf.command.out
    cp .command.run updater/logs/nf.command.run
    cp .command.sh updater/logs/nf.command.sh
    cp .command.trace updater/logs/nf.command.trace

    cat <<-END_VERSIONS > updater/logs/versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
        amrfinderplus-database: \$(echo \$(echo \$(amrfinder --database amrfinderplus --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev))
    END_VERSIONS
    """
}
