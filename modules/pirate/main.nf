/**
 * Pangenome Identification and Reconciliation Analysis Tool for Epidemiology (PIRATE).
 *
 * Uses [PIRATE](https://github.com/SionBayliss/PIRATE) to construct the pangenome of a
 * collection of bacterial isolates. It clusters orthologous genes and generates the core
 * genome alignment and a gene presence/absence matrix, which is compatible with downstream
 * analysis tools like Scoary for association testing.
 *
 * @status stable
 * @keywords pan-genome, orthologs, core genome, gene presence-absence, epidemiology, annotation
 * @tags complexity:complex input-type:multiple output-type:multiple features:compression
 * @citation pirate
 *
 * @input record(meta, gff)
 * - `meta`: Groovy Map containing sample information
 * - `gff`: A list of annotated genome files in GFF3 format
 *
 * @output record(meta, aln, csv, results, logs, nf_logs, versions)
 * - `aln`: The core-genome alignment (*core-genome.aln.gz), suitable for phylogenetic tree building
 * - `csv`: Gene presence/absence matrix in CSV format, compatible with Scoary
 *
 * @results pirate/
 * - `pirate/*`: Full PIRATE output directory (co-ords, modified GFFs, allele files, etc.)
 */
nextflow.preview.types = true

process PIRATE {
    tag "${prefix}"
    label 'process_high'
    label 'process_long'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, gff: Set<Path>): Record

    stage:
    stageAs 'gff-tmp/*', gff

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        aln: file("${prefix}.aln.gz", optional: true),
        csv: file("pirate/gene_presence_absence.csv", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.aln.gz", optional: true),
            files("pirate/gene_presence_absence.csv", optional: true),
            files("pirate/*")
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

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
        cp supplemental/core_alignment.fasta.gz ./${prefix}.aln.gz
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
