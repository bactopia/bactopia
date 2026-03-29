/**
 * Fast and scalable bacterial pangenome analysis using a graph-based approach.
 *
 * Uses [Panaroo](https://gtonkinhill.github.io/panaroo/) to cluster genes from multiple
 * annotated bacterial genomes into orthologous groups, correcting for gene splitting and
 * merges. The primary outputs are the gene presence/absence matrix (the pan-genome) and
 * a core-genome alignment (for phylogenetics).
 *
 * @status stable
 * @keywords pan-genome, orthologs, core genome, gene presence-absence, graph-based, annotation
 * @tags complexity:complex input-type:multiple output-type:multiple features:compression
 * @citation panaroo
 *
 * @input record(meta, gff)
 * - `meta`: Groovy Map containing sample information
 * - `gff`: A list of annotated genome files in GFF3 format (required input)
 *
 * @output record(meta, aln, filtered_aln, csv, panaroo_csv, results, logs, nf_logs, versions)
 * - `aln`: The core-genome alignment (*core-genome.aln.gz), suitable for phylogenetic tree building
 * - `filtered_aln`: The core-genome alignment with highly recombinant regions filtered out (optional)
 * - `csv`: Gene presence/absence matrix in Roary-compatible CSV format (optional)
 * - `panaroo_csv`: Gene presence/absence matrix in Panaroo's native CSV format (optional)
 *
 * @results panaroo/
 * - `panaroo/*`: Full Panaroo output directory (graphs, gene data, struct files, etc.)
 */
nextflow.preview.types = true

process PANAROO_RUN {
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
        filtered_aln: file("${prefix}.filtered.aln.gz", optional: true),
        csv: file("panaroo/gene_presence_absence_roary.csv", optional: true),
        panaroo_csv: file("panaroo/gene_presence_absence.csv", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.aln.gz", optional: true),
            files("${prefix}.filtered.aln.gz", optional: true),
            files("panaroo/gene_presence_absence_roary.csv", optional: true),
            files("panaroo/gene_presence_absence.csv", optional: true),
            files("panaroo/*")
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
    meta.output_dir = ""
    meta.logs_dir = "panaroo/logs"
    meta.process_name = task.ext.process_name
    """
    mkdir gff
    cp -L gff-tmp/* gff/
    find gff/ -name "*.gz" | xargs gunzip

    # Make FOFN of gff (Prokka) and gff3 (Bakta) files
    find gff/ -name "*.gff" -or -name "*.gff3" > gff-fofn.txt

    panaroo \\
        ${task.ext.args} \\
        -t ${task.cpus} \\
        -o supplemental \\
        -i gff-fofn.txt

    # Cleanup
    find . -name "*.fas" | xargs -I {} -P ${task.cpus} -n 1 gzip {}
    find . -name "*.fa" | xargs -I {} -P ${task.cpus} -n 1 gzip {}
    find . -name "*.fasta" | xargs -I {} -P ${task.cpus} -n 1 gzip {}
    find . -name "*.aln" | xargs -I {} -P ${task.cpus} -n 1 gzip {}
    find . -name "*.gml" | xargs -I {} -P ${task.cpus} -n 1 gzip {}

    if [[ -f "supplemental/core_gene_alignment.aln.gz" ]]; then
        mv supplemental/core_gene_alignment.aln.gz ./${prefix}.aln.gz
    fi

    if [[ -f "supplemental/core_gene_alignment_filtered.aln.gz" ]]; then
        mv supplemental/core_gene_alignment_filtered.aln.gz ./${prefix}.filtered.aln.gz
    fi

    if [[ -f "supplemental/gene_data.csv" ]]; then
        gzip supplemental/gene_data.csv
    fi

    mv supplemental/ panaroo/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panaroo: \$(echo \$(panaroo --version 2>&1) | sed 's/^.*panaroo //')
    END_VERSIONS
    """
}
