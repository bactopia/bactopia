process MCRONI {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv")         , emit: tsv
    tuple val(meta), path("*.fa")          , emit: fa, optional: true
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
    def VERSION = '1.0.4' // Version information not provided by tool on CLI
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    mcroni \\
        --output $prefix \\
        --fasta $fasta_name

    EX_COLS=\$(head -n 1 ${prefix}_table.tsv| tr '\\t' '\\n' | wc -l)
    OBS_COLS=\$(tail -n 1 ${prefix}_table.tsv| tr '\\t' '\\n' | wc -l)
    if [ "\$EX_COLS" != "\$OBS_COLS" ]; then
        sed -i 's/NA\$/NA\\tNA/' ${prefix}_table.tsv
    fi

    # Cleanup
    rm -rf ${fasta_name} ${fasta_name}.ndb ${fasta_name}.not ${fasta_name}.ntf ${fasta_name}.nto

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mcroni: $VERSION
    END_VERSIONS
    """
}
