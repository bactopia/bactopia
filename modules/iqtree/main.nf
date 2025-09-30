process IQTREE {
    tag "${prefix}"
    label 'process_medium'
    label 'process_long'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(alignment)

    output:
    tuple val(meta), path("${process_name}/*")      , emit: supplemental
    tuple val(meta), path(treefile)                 , emit: phylogeny
    tuple val(meta), path(alignment)                , emit: alignment
    tuple val(meta), path(alignment), path(treefile), emit: aln_tree
    tuple val(meta), path("*.{log,err}")            , emit: logs, optional: true
    tuple val(meta), path(".command.begin") , emit: nf_begin
    tuple val(meta), path(".command.err")   , emit: nf_err
    tuple val(meta), path(".command.log")   , emit: nf_log
    tuple val(meta), path(".command.out")   , emit: nf_out
    tuple val(meta), path(".command.run")   , emit: nf_run
    tuple val(meta), path(".command.sh")    , emit: nf_sh
    tuple val(meta), path(".command.trace") , emit: nf_trace
    tuple val(meta), path("versions.yml")   , emit: versions

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"
    process_name = _meta.process_name == "iqtree-fast" ? "iqtree-fast" : task.ext.process_name
    args = process_name == "iqtree-fast" ? task.ext.fast_args : task.ext.args
    treefile = process_name == "iqtree-fast" ? "${process_name}/${prefix}.treefile" : "${prefix}.treefile"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.output_dir = "${task.ext.rundir}/"
    meta.logs_dir = "${task.ext.rundir}/${process_name}/logs/"
    meta.process_name = process_name
    """
    iqtree \\
        $args \\
        -s $alignment \\
        -nt $task.cpus \\
        -ntmax $task.cpus \\
        -pre $prefix

    # Only gzip files if they exist
    if [[ -f "${prefix}.alninfo" ]]; then
        gzip ${prefix}.alninfo
    fi

    mkdir temp
    mv ${prefix}* temp/
    mv temp/ ${process_name}/

    if [ "${process_name}" != "iqtree-fast" ]; then
        mv ${process_name}/${prefix}.treefile ./
        mv ${process_name}/${alignment} ./
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
    END_VERSIONS
    """
}
