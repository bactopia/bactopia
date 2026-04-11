/**
 * Identify genes by local assembly of reads.
 *
 * Uses [ARIBA](https://github.com/sanger-pathogens/ariba) (Antimicrobial Resistance Identification
 * By Assembly) to detect AMR and virulence genes by creating local assemblies from paired-end reads.
 *
 * Uses explicit positional record fields for reads:
 * - Input: record(meta, r1, r2) where each read slot is Path
 *
 * @status stable
 * @keywords fastq, local assembly, antimicrobial resistance, virulence, ariba
 * @tags complexity:moderate input-type:multiple output-type:multiple features:archive-output,compression
 * @citation ariba
 *
 * @input record(meta, r1, r2)
 * - `meta`: Groovy Record containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 *
 * @input db
 * An [ARIBA](https://github.com/sanger-pathogens/ariba) prepared database
 *
 * @output record(meta, report, summary, results, logs, nf_logs, versions)
 * - `report`: Tab-delimited detailed report of gene detection results
 * - `summary`: Comma-separated condensed summary of detected genes
 *
 * @results supplemental
 * - `assembled_genes.fa.gz`: FASTA of locally assembled gene sequences
 * - `assembled_seqs.fa.gz`: FASTA of assembled contig sequences
 * - `assemblies.fa.gz`: Raw assembly output sequences
 * - `log.clusters.gz`: Per-cluster log of assembly and mapping details
 * - `version_info.txt`: ARIBA version and database information
 */
nextflow.preview.types = true

process ARIBA_RUN {
    tag "${prefix} - ${db.name}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Record,
        r1: Path,
        r2: Path
    )
    db: Path

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        report: file("${prefix}-report.tsv"),
        summary: file("${prefix}-summary.csv"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}-report.tsv"),
            files("${prefix}-summary.csv"),
            files("supplemental/*")
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
        output_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}",
        logs_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}",
        process_name: task.ext.process_name
    )

    def db_name = db.getName().replace('.tar.gz', '')
    """
    tar -xzvf ${db}
    mv ${db_name} ${db_name}db
    ariba \\
        run \\
        ${db_name}db/ \\
        ${r1} ${r2} \\
        ${db_name} \\
        ${task.ext.args} \\
        --threads ${task.cpus}

    ariba \\
        summary \\
        ${db_name}/summary \\
        ${db_name}/report.tsv \\
        --cluster_cols assembled,match,known_var,pct_id,ctg_cov,novel_var \\
        --col_filter n \\
        --row_filter n

    # Rename to avoid naming collisions
    mv ${db_name}/report.tsv ./${prefix}-report.tsv
    mv ${db_name}/summary.csv ./${prefix}-summary.csv
    mv ${db_name}/ supplemental/

    # Cleanup
    rm -rf ${db_name}db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ariba:  \$(echo \$(ariba version 2>&1) | sed 's/^.*ARIBA version: //;s/ .*\$//')
    END_VERSIONS
    """
}
