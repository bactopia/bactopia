/**
 * Pangenome toolbox for bacterial genomes.
 *
 * This process executes pirate to perform analysis
 *
 * @status stable
 * @keywords gff, pan-genome, alignment
 * @tags complexity:moderate input-type:single output-type:multiple features:archive-output, compression, conditional-logic
 * @citation pirate
 *
 * @input tuple(meta, gff)
 * - `meta`: Groovy Map containing sample information
 * - `gff`: A set of GFF3 formatted files
 *
 * @output supplemental Supplemental
 * @output aln          Core-genome alignment produced by PIRATE (Optional)
 * @output csv          Gene presence/absence CSV compatible with Scoary
 * @output logs         Optional tool execution logs
 * @output nf_logs      Nextflow execution logs
 * @output versions     Software version information (YAML format)
 */
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
    supplemental = tuple(meta, files("pirate/*"))
    aln          = tuple(meta, file("core-genome.aln.gz", optional: true))
    csv          = tuple(meta, file("pirate/gene_presence_absence.csv", optional: true))
    logs         = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs      = tuple(meta, files(".command.*"))
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
