process MASHTREE {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(seqs)

    output:
    tuple val(meta), path("*.dnd")         , emit: tree
    tuple val(meta), path("*.tsv")         , emit: matrix
    tuple val(meta), path("sketches/*")    , emit: sketches, optional: true
    tuple val(meta), path("*.{log,err}")   , emit: logs, optional: true
    tuple val(meta), path(".command.begin"), emit: nf_begin
    tuple val(meta), path(".command.err")  , emit: nf_err
    tuple val(meta), path(".command.log")  , emit: nf_log
    tuple val(meta), path(".command.out")  , emit: nf_out
    tuple val(meta), path(".command.run")  , emit: nf_run
    tuple val(meta), path(".command.sh")   , emit: nf_sh
    tuple val(meta), path(".command.trace"), emit: nf_trace
    tuple val(meta), path("versions.yml")  , emit: versions

    script:
    prefix = task.ext.prefix ?: "${_meta.id}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = ""
    meta.logs_dir = "logs/"
    meta.process_name = task.ext.process_name
    """
    mkdir mashtree-tmp

    mashtree \\
        ${task.ext.args} \\
        --numcpus $task.cpus \\
        --outmatrix ${prefix}.tsv \\
        --outtree ${prefix}.dnd \\
        --tempdir mashtree-tmp/ \\
        $seqs

    # Clean up
    rm -rf mashtree-tmp/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(echo \$(mash 2>&1) | sed 's/^.*Mash version //;s/ .*\$//')
        mashtree: \$( echo \$( mashtree --version 2>&1 ) | sed 's/^.*Mashtree //' )
    END_VERSIONS
    """
}
