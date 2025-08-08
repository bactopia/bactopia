process MCRONI {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}" :
        "${task.ext.docker}" }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv")               , emit: tsv
    tuple val(meta), path("*.fa"), optional: true, emit: fa
    path "versions.yml"                           , emit: versions
    path ".command.begin"                         , emit: begin
    path ".command.err"                           , emit: err
    path ".command.log"                           , emit: log
    path ".command.out"                           , emit: out
    path ".command.run"                           , emit: run
    path ".command.sh"                            , emit: sh
    path ".command.trace"                         , emit: trace

    script:
    def VERSION = '1.0.4' // Version information not provided by tool on CLI
    prefix = task.ext.prefix ?: "${meta.id}"
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
