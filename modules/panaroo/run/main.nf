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
 * - `meta`: Groovy Record containing sample information
 * - `gff`: A list of annotated genome files in GFF3 format (required input)
 *
 * @output record(meta, aln?, filtered_aln?, csv?, panaroo_csv?, results, logs, nf_logs, versions)
 * - `aln?`: The core-genome alignment (*core-genome.aln.gz), suitable for phylogenetic tree building
 * - `filtered_aln?`: The core-genome alignment with highly recombinant regions filtered out
 * - `csv?`: Gene presence/absence matrix in Roary-compatible CSV format
 * - `panaroo_csv?`: Gene presence/absence matrix in Panaroo's native CSV format
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
    record (
        meta: Record,
        gff: Set<Path>
    )

    stage:
    stageAs gff, 'staging/gff/*'

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
    def _meta = meta
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = record(
        id: "${prefix}-${task.process}",
        name: prefix,
        scope: task.ext.scope,
        output_dir: "",
        logs_dir: "panaroo/logs",
        process_name: task.ext.process_name
    )
    """
    mkdir gff
    cp -L staging/gff/* gff/
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
    rm -rf gff/ gff-fofn.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panaroo: \$(echo \$(panaroo --version 2>&1) | sed 's/^.*panaroo //')
    END_VERSIONS
    """
}
