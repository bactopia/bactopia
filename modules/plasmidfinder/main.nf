process PLASMIDFINDER {
    tag "$meta.id"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.json")                    , emit: json
    tuple val(meta), path("*.txt")                     , emit: txt
    tuple val(meta), path("${prefix}.tsv")             , emit: tsv
    tuple val(meta), path("*-hit_in_genome_seq.fsa.gz"), emit: genome_seq
    tuple val(meta), path("*-plasmid_seqs.fsa.gz")     , emit: plasmid_seq
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
    def args = task.ext.args ?: ''
    def VERSION = '2.1.6' // Version information not provided by tool on CLI
    prefix = task.ext.prefix ?: "${meta.id}"
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
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
