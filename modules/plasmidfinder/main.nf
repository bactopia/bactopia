process PLASMIDFINDER {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}" :
        "${task.ext.docker}" }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.json")                   , emit: json
    tuple val(meta), path("*.txt")                    , emit: txt
    tuple val(meta), path("${prefix}.tsv")            , emit: tsv
    tuple val(meta), path("*-hit_in_genome_seq.fsa.gz"), emit: genome_seq
    tuple val(meta), path("*-plasmid_seqs.fsa.gz")     , emit: plasmid_seq
    path "versions.yml"                               , emit: versions
    path ".command.begin"                             , emit: begin
    path ".command.err"                               , emit: err
    path ".command.log"                               , emit: log
    path ".command.out"                               , emit: out
    path ".command.run"                               , emit: run
    path ".command.sh"                                , emit: sh
    path ".command.trace"                             , emit: trace

    script:
    def args = task.ext.args ?: ''
    def VERSION = '2.1.6' // Version information not provided by tool on CLI
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    plasmidfinder.py \\
        $args \\
        -i $fasta_name \\
        -o ./ \\
        -x

    # Rename hard-coded outputs with prefix to avoid name collisions
    mv data.json ${prefix}.json
    mv results.txt ${prefix}.txt
    mv Hit_in_genome_seq.fsa ${prefix}-hit_in_genome_seq.fsa
    mv Plasmid_seqs.fsa ${prefix}-plasmid_seqs.fsa

    # Add sample name to TSV results
    head -n 1 results_tab.tsv | sed "s/^/Sample\t/" > ${prefix}.tsv
    tail -n +2 results_tab.tsv | sed "s/^/${prefix}\t/" >> ${prefix}.tsv

    # Cleanup
    gzip *.fsa
    rm -rf ${fasta_name} results_tab.tsv tmp/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plasmidfinder: $VERSION
    END_VERSIONS
    """
}
