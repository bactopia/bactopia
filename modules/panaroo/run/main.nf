process PANAROO_RUN {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'

    conda "${task.ext.conda}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}" :
        "${task.ext.docker}" }"

    input:
    tuple val(meta), path(gff, stageAs: 'gff-tmp/*')

    output:
    tuple val(meta), path("results/*")                                              , emit: results
    tuple val(meta), path("core-genome.aln.gz")                     , optional: true, emit: aln
    tuple val(meta), path("results/gene_presence_absence_roary.csv"), optional: true, emit: csv
    tuple val(meta), path("results/gene_presence_absence.csv")      , optional: true, emit: panaroo_csv
    path "versions.yml"                                                             , emit: versions
    path ".command.begin"                                                           , emit: begin
    path ".command.err"                                                             , emit: err
    path ".command.log"                                                             , emit: log
    path ".command.out"                                                             , emit: out
    path ".command.run"                                                             , emit: run
    path ".command.sh"                                                              , emit: sh
    path ".command.trace"                                                           , emit: trace

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir gff
    cp -L gff-tmp/* gff/
    find gff/ -name "*.gz" | xargs gunzip

    # Make FOFN of gff (Prokka) and gff3 (Bakta) files
    find gff/ -name "*.gff" -or -name "*.gff3" > gff-fofn.txt

    panaroo \\
        $args \\
        -t $task.cpus \\
        -o results \\
        -i gff-fofn.txt

    # Cleanup
    find . -name "*.fas" | xargs -I {} -P $task.cpus -n 1 gzip {}
    find . -name "*.fa" | xargs -I {} -P $task.cpus -n 1 gzip {}
    find . -name "*.fasta" | xargs -I {} -P $task.cpus -n 1 gzip {}
    find . -name "*.aln" | xargs -I {} -P $task.cpus -n 1 gzip {}
    find . -name "*.gml" | xargs -I {} -P $task.cpus -n 1 gzip {}

    if [[ -f "results/core_gene_alignment.aln.gz" ]]; then
        cp results/core_gene_alignment.aln.gz ./core-genome.aln.gz
    fi

    if [[ -f "results/gene_data.csv" ]]; then
        gzip results/gene_data.csv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panaroo: \$(echo \$(panaroo --version 2>&1) | sed 's/^.*panaroo //')
    END_VERSIONS
    """
}
