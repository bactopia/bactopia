process GAMMA {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    tuple val(meta), path("*.gamma")               , emit: gamma
    tuple val(meta), path("*.psl")                 , emit: psl
    tuple val(meta), path("*.gff")  , optional:true, emit: gff
    tuple val(meta), path("*.fasta"), optional:true, emit: fasta
    tuple val(meta), path("*.{log,err}")           , emit: logs, optional: true
    tuple val(meta), path(".command.begin")        , emit: nf_begin
    tuple val(meta), path(".command.err")          , emit: nf_err
    tuple val(meta), path(".command.log")          , emit: nf_log
    tuple val(meta), path(".command.out")          , emit: nf_out
    tuple val(meta), path(".command.run")          , emit: nf_run
    tuple val(meta), path(".command.sh")           , emit: nf_sh
    tuple val(meta), path(".command.trace")        , emit: nf_trace
    tuple val(meta), path("versions.yml")          , emit: versions

    script:
    prefix = task.ext.prefix ? "${task.ext.prefix}" : "${meta.id}"
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    def VERSION = '2.1' // Version information not provided by tool on CLI
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    GAMMA.py \\
        $task.ext.args \\
        $fasta_name \\
        $db \\
        $prefix

    # Cleanup
    rm -rf ${fasta_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gamma: $VERSION
    END_VERSIONS
    """
}
