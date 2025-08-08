process MOBSUITE_RECON {
    tag "$meta.id"
    label 'process_medium'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}" :
        "${task.ext.docker}" }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("results/chromosome.fasta.gz")    , emit: chromosome
    tuple val(meta), path("results/contig_report.txt")      , emit: contig_report
    tuple val(meta), path("results/plasmid_*.fasta.gz")     , emit: plasmids        , optional: true
    tuple val(meta), path("results/${prefix}-mobtyper.txt") , emit: mobtyper_results, optional: true
    path "versions.yml"                                      , emit: versions
    path ".command.begin"                                    , emit: begin
    path ".command.err"                                      , emit: err
    path ".command.log"                                      , emit: log
    path ".command.out"                                      , emit: out
    path ".command.run"                                      , emit: run
    path ".command.sh"                                       , emit: sh
    path ".command.trace"                                    , emit: trace

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    mob_recon \\
        --infile $fasta_name \\
        $args \\
        --num_threads $task.cpus \\
        --outdir results \\
        --sample_id $prefix

    if [[ -f "results/mobtyper_results.txt" ]]; then
        mv results/mobtyper_results.txt results/${prefix}-mobtyper.txt
    fi

    # Cleanup
    gzip results/*.fasta
    rm -rf ${fasta_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mobsuite: \$(echo \$(mob_recon --version 2>&1) | sed 's/^.*mob_recon //; s/ .*\$//')
    END_VERSIONS
    """
}
