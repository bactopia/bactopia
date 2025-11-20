nextflow.preview.types = true

process PIRATE {
    tag "${prefix}"
    label 'process_high'
    label 'process_long'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, gff) : Tuple<Map, Set<Path>>

    stage:
    stageAs 'gff-tmp/*', gff

    output:
    supplemental = tuple(meta, file("pirate/*"))
    aln          = tuple(meta, file("core-genome.aln.gz", optional: true))
    csv          = tuple(meta, file("pirate/gene_presence_absence.csv", optional: true))
    logs         = tuple(meta, file("*.{log,err}", optional: true))
    nf_begin     = tuple(meta, file(".command.begin"))
    nf_err       = tuple(meta, file(".command.err"))
    nf_log       = tuple(meta, file(".command.log"))
    nf_out       = tuple(meta, file(".command.out"))
    nf_run       = tuple(meta, file(".command.run"))
    nf_sh        = tuple(meta, file(".command.sh"))
    nf_trace     = tuple(meta, file(".command.trace"))
    versions     = tuple(meta, file("versions.yml"))

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.process_name = task.ext.process_name
    meta.output_dir = ""
    meta.logs_dir = "${meta.process_name}/logs"
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
        --threads ${task.cpus} \\
        --input ./gff/ \\
        --output supplemental/
    PIRATE_to_roary.pl -i supplemental/PIRATE.*.tsv -o supplemental/gene_presence_absence.csv
    find . -name "*.fasta" | xargs -I {} -P ${task.cpus} -n 1 gzip {}

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
