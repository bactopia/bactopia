process PIRATE {
    tag "${prefix}"
    label 'process_high'
    label 'process_long'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(gff, stageAs: 'gff-tmp/*')

    output:
    tuple val(meta), path("pirate/*")                        , emit: results
    tuple val(meta), path("core-genome.aln.gz")              , emit: aln, optional: true
    tuple val(meta), path("pirate/gene_presence_absence.csv"), emit: csv, optional: true
    tuple val(meta), path("pirate/gene_presence_absence.csv"), emit: panaroo_csv, optional: true
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
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.output_dir = "${task.ext.rundir}/"
    meta.logs_dir = "${task.ext.rundir}/${task.ext.process_name}/logs"
    meta.process_name = task.ext.process_name
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
        ${task.ext.args} \\
        --align \\
        --threads $task.cpus \\
        --input ./gff/ \\
        --output supplemental/
    PIRATE_to_roary.pl -i supplemental/PIRATE.*.tsv -o supplemental/gene_presence_absence.csv
    find . -name "*.fasta" | xargs -I {} -P $task.cpus -n 1 gzip {}

    # Only copy files if they exist
    if [[ -f "supplemental/core_alignment.fasta.gz" ]]; then
        cp supplemental/core_alignment.fasta.gz ./core-genome.aln.gz
    fi

    # Cleanup
    rm -rf gff/
    gzip supplemental/co-ords/*.tab
    gzip supplemental/modified_gffs/*.gff

    mv supplemental/ pirate/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pirate: \$( echo \$( PIRATE --version 2>&1) | sed 's/PIRATE //' )
    END_VERSIONS
    """
}
