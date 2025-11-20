nextflow.preview.types = true

process BAKTA_DOWNLOAD {
    label 'process_low'
    label 'process_long'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    output:
    db         = file("bakta-${task.ext.bakta_db_type}/*", optional: true)
    db_tarball = file("bakta-${task.ext.bakta_db_type}.tar.gz", optional: true)
    logs       = file("logs/*", optional: true)

    script:
    """
    bakta_db \\
        download \\
        --type "${task.ext.bakta_db_type}" \\
        --output bakta

    if [ "${task.ext.save_as_tarball}" == "true" ]; then
        tar -czf bakta-${task.ext.bakta_db_type}.tar.gz bakta/
        rm -rf bakta/
    else
        mv bakta/ bakta-${task.ext.bakta_db_type}/
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
        bakta: \$( echo \$(bakta --version 2>&1) | sed 's/^.*bakta //' )
    END_VERSIONS
    """
}
