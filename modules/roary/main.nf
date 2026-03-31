/**
 * Rapid large-scale prokaryote pan genome analysis.
 *
 * Uses [Roary](https://github.com/sanger-pathogens/Roary) to calculate the pan genome of a
 * collection of prokaryotic annotated assemblies. It outputs a core gene alignment and a
 * gene presence/absence table.
 *
 * @status stable
 * @keywords pangenome, orthology, core genome, alignment, bacteria
 * @tags complexity:moderate input-type:multiple output-type:multiple features:conditional-logic
 * @citation roary
 *
 * @input record(meta, gff)
 * - `meta`: Groovy Map containing sample information
 * - `gff`: List of GFF3 files to be analyzed (typically from Prokka)
 *
 * @output record(meta, aln, csv, results, logs, nf_logs, versions)
 * - `aln`: Core genome alignment in FASTA format
 * - `csv`: Gene presence/absence table
 *
 * @results roary/
 * - `roary/*`: Full Roary output directory (cluster files, alignments, Rtab, etc.)
 */
nextflow.preview.types = true

process ROARY {
    tag "${prefix}"
    label 'process_high'
    label 'process_long'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (meta: Map, gff: Set<Path>): Record

    stage:
    stageAs 'gff-tmp/*', gff

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        aln: file("${prefix}.aln.gz", optional: true),
        csv: file("roary/gene_presence_absence.csv", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.aln.gz", optional: true),
            files("roary/gene_presence_absence.csv", optional: true),
            files("roary/*")
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    def _meta = meta
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

    # Roary only supports .gff extension, will need to adjust for gff3 (Bakta) files
    # https://github.com/sanger-pathogens/Roary/blob/master/lib/Bio/Roary/PrepareInputFiles.pm#L82
    # note for later: "xargs -r" will not run if no files are found
    find gff/ -name "*.gff3.gz" | xargs -r gunzip
    find gff/ -name "*.gff3" -print0 | while read -d \$'\0' file; do mv "\$file" "\${file%.gff3}.gff"; done

    roary \\
        ${task.ext.args} \\
        -p ${task.cpus} \\
        -f supplemental/ \\
        gff/*.gff

    gzip supplemental/*.aln
    gzip supplemental/*.fa

    if [[ -f "supplemental/core_gene_alignment.aln.gz" ]]; then
        cp supplemental/core_gene_alignment.aln.gz ./${prefix}.aln.gz
    fi

    # Cleanup
    rm -rf gff/
    mv supplemental/ roary/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        roary: \$( roary --version )
    END_VERSIONS
    """
}
