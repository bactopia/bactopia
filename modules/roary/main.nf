process ROARY {
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
    path "versions.yml", emit: versions
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
    """
    mkdir gff
    cp -L gff-tmp/* gff/
    find gff/ -name "*.gff.gz" | xargs -r gunzip

    # Roary only supports .gff extension, will need to adjust for gff3 (Bakta) files
    # https://github.com/sanger-pathogens/Roary/blob/master/lib/Bio/Roary/PrepareInputFiles.pm#L82
    # note for later: "xargs -r" will not run if no files are found
    find gff/ -name "*.gff3.gz" | xargs -r gunzip
    find gff/ -name "*.gff3" -print0 | while read -d \$'\0' file; do mv "\$file" "\${file%.gff3}.gff"; done

    roary \\
        $args \\
        -p $task.cpus \\
        -f results/ \\
        gff/*.gff

    gzip results/*.aln
    gzip results/*.fa

    if [[ -f "results/core_gene_alignment.aln.gz" ]]; then
        cp results/core_gene_alignment.aln.gz ./core-genome.aln.gz
    fi

    # clean up
    rm -rf gff/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        roary: \$( roary --version )
    END_VERSIONS
    """
}
