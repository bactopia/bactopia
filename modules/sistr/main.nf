process SISTR {
    tag "$meta.id"
    label 'process_medium'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}" :
        "${task.ext.docker}" }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv")            , emit: tsv
    tuple val(meta), path("*-allele.fasta.gz"), emit: allele_fasta
    tuple val(meta), path("*-allele.json.gz") , emit: allele_json
    tuple val(meta), path("*-cgmlst.csv")     , emit: cgmlst_csv
    path "versions.yml"                       , emit: versions
    path ".command.begin", emit: begin
    path ".command.err", emit: err
    path ".command.log", emit: log
    path ".command.out", emit: out
    path ".command.run", emit: run
    path ".command.sh", emit: sh
    path ".command.trace", emit: trace

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    sistr \\
        --qc \\
        $args \\
        --threads $task.cpus \\
        --alleles-output ${prefix}-allele.json \\
        --novel-alleles ${prefix}-allele.fasta \\
        --cgmlst-profiles ${prefix}-cgmlst.csv \\
        --output-prediction ${prefix} \\
        --output-format tab \\
        $fasta_name

    mv ${prefix}.tab ${prefix}.tsv
    gzip ${prefix}-allele.json
    gzip ${prefix}-allele.fasta

    # Cleanup
    rm -rf ${fasta_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sistr: \$(echo \$(sistr --version 2>&1) | sed 's/^.*sistr_cmd //; s/ .*\$//' )
    END_VERSIONS
    """
}
