process PANAROO_RUN {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(meta), path(gff, stageAs: 'gff-tmp/*')

    output:
    tuple val(meta), path("results/*")                                              , emit: results
    tuple val(meta), path("core-genome.aln.gz")                     , optional: true, emit: aln
    tuple val(meta), path("results/gene_presence_absence_roary.csv"), optional: true, emit: csv
    tuple val(meta), path("results/gene_presence_absence.csv")      , optional: true, emit: panaroo_csv
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
    prefix = task.ext.prefix ?: "${meta.id}"
    meta.output_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${meta.id}/tools/${task.ext.process_name}/${task.ext.subdir}/logs"
    meta.process_name = task.ext.process_name
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
