process AGRVATE {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(fasta, stageAs: 'input/*')

    output:
    tuple val(meta), path("${prefix}-summary.tab"), emit: summary
    tuple val(meta), path("supplemental/*")       , emit: supplemental
    tuple val(meta), path("*.{log,err}")          , emit: logs, optional: true
    tuple val(meta), path(".command.begin")       , emit: nf_begin
    tuple val(meta), path(".command.err")         , emit: nf_err
    tuple val(meta), path(".command.log")         , emit: nf_log
    tuple val(meta), path(".command.out")         , emit: nf_out
    tuple val(meta), path(".command.run")         , emit: nf_run
    tuple val(meta), path(".command.sh")          , emit: nf_sh
    tuple val(meta), path(".command.trace")       , emit: nf_trace
    tuple val(meta), path("versions.yml")         , emit: versions

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = "${prefix}.fna"
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > ./$fasta_name
    else
        cat $fasta > ./$fasta_name
    fi

    agrvate \\
        $task.ext.args \\
        -i $fasta_name

    mv ${meta.name}-results/ supplemental/
    mv supplemental/${meta.name}-summary.tab ./

    # Cleanup
    rm -rf ${fasta_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agrvate: \$(echo \$(agrvate -v 2>&1) | sed 's/agrvate v//;')
    END_VERSIONS
    """
}
