process STAPHOPIASCCMEC {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.conda_env}"
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/staphopia-sccmec:1.0.0--hdfd78af_0' :
        'quay.io/biocontainers/staphopia-sccmec:1.0.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "*.{log,err}", emit: logs, optional: true
    path ".command.begin"   , emit: begin
    path ".command.err"     , emit: err
    path ".command.log"     , emit: log
    path ".command.out"     , emit: out
    path ".command.run"     , emit: run
    path ".command.sh"      , emit: sh
    path ".command.trace"   , emit: trace
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    staphopia-sccmec --assembly ${fasta_name} ${task.ext.args} > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        staphopia-sccmec: \$(staphopia-sccmec --version 2>&1 | sed 's/^.*staphopia-sccmec //')
    END_VERSIONS
    """
}
