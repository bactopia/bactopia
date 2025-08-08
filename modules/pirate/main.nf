process PIRATE {
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
    tuple val(meta), path("results/*")                        , emit: results
    tuple val(meta), path("core-genome.aln.gz")               , emit: aln, optional: true
    tuple val(meta), path("results/gene_presence_absence.csv"), emit: csv, optional: true
    path "versions.yml"                                       , emit: versions
    path ".command.begin"                                     , emit: begin
    path ".command.err"                                       , emit: err
    path ".command.log"                                       , emit: log
    path ".command.out"                                       , emit: out
    path ".command.run"                                       , emit: run
    path ".command.sh"                                        , emit: sh
    path ".command.trace"                                     , emit: trace

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir gff
    cp -L gff-tmp/* gff/
    find gff/ -name "*.gff.gz" | xargs -r gunzip

    # PIRATE only supports .gff extension, will need to adjust for gff3 (Bakta) files
    # https://github.com/SionBayliss/PIRATE/blob/master/scripts/run_PIRATE.pl#L153
    # note for later: "xargs -r" will not run if no files are found
    find gff/ -name "*.gff3.gz" | xargs -r gunzip
    find gff/ -name "*.gff3" -print0 | while read -d \$'\0' file; do mv "\$file" "\${file%.gff3}.gff"; done

    PIRATE \\
        $args \\
        --align \\
        --threads $task.cpus \\
        --input ./gff/ \\
        --output results/
    PIRATE_to_roary.pl -i results/PIRATE.*.tsv -o results/gene_presence_absence.csv
    find . -name "*.fasta" | xargs -I {} -P $task.cpus -n 1 gzip {}

    # Only copy files if they exist
    if [[ -f "results/core_alignment.fasta.gz" ]]; then
        cp results/core_alignment.fasta.gz ./core-genome.aln.gz
    fi

    # Cleanup
    rm -rf gff/
    gzip results/co-ords/*.tab
    gzip results/modified_gffs/*.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pirate: \$( echo \$( PIRATE --version 2>&1) | sed 's/PIRATE //' )
    END_VERSIONS
    """
}
