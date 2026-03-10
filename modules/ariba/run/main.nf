/**
 * Identify genes by local assembly of reads.
 *
 * Uses [ARIBA](https://github.com/sanger-pathogens/ariba) (Antimicrobial Resistance Identification
 * By Assembly) to detect AMR and virulence genes by creating local assemblies from paired-end reads.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords fastq, local assembly, antimicrobial resistance, virulence, ariba
 * @tags complexity:moderate input-type:multiple output-type:multiple features:archive-output,compression
 * @citation ariba
 *
 * @input record(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads (not supported by ARIBA)
 * - `lr`: Long reads (not supported by ARIBA)
 *
 * @input db
 * An [ARIBA](https://github.com/sanger-pathogens/ariba) prepared database
 *
 * @output record(meta, report, summary, results, logs, nf_logs, versions)
 * - `report`: Tab-delimited detailed report of gene detection results
 * - `summary`: Comma-separated condensed summary of detected genes
 */
nextflow.preview.types = true

process ARIBA_RUN {
    tag "${_meta.name} - ${db.name}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, r1: Path?, r2: Path?, se: Path?, lr: Path?): Record
    db: Path

    output:
    record(
        meta: meta,
        // Named fields (upstream consumers access these)
        report: file("${prefix}-report.tsv"),
        summary: file("${prefix}-summary.csv"),
        // Generic fields (same convention across every module)
        results: files("${prefix}-report.tsv") + files("${prefix}-summary.csv") + files("supplemental/*"),
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
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name

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
