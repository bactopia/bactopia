nextflow.preview.types = true

process WGET {
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    _meta : Map

    output:
    download = file("${prefix}/${filename}")
    logs     = files("${prefix}/logs/*", optional: true)

    script:
    prefix = _meta.name
    filename = _meta.save_as
    """
    wget ${task.ext.args} -O ${filename} ${_meta.url}

    # Move outputs to tool specific folder
    mkdir -p ${prefix}/logs
    mv ${filename} ${prefix}/
    cp .command.begin ${prefix}/logs/nf.command.begin
    cp .command.err ${prefix}/logs/nf.command.err
    cp .command.log ${prefix}/logs/nf.command.log
    cp .command.out ${prefix}/logs/nf.command.out
    cp .command.run ${prefix}/logs/nf.command.run
    cp .command.sh ${prefix}/logs/nf.command.sh
    cp .command.trace ${prefix}/logs/nf.command.trace

    cat <<-END_VERSIONS > ${prefix}/logs/versions.yml
    "${task.process}":
        wget: \$(echo \$( wget --version 2>&1) | sed "s/GNU Wget //;s/ .*//;" )
    END_VERSIONS
    """
}
