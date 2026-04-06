/**
 * k-mer based Streptococcus pneumoniae serotyping.
 *
 * Uses [SeroBA](https://github.com/sanger-pathogens/seroba) to identify the serotype of
 * *Streptococcus pneumoniae* from Illumina paired-end reads using a k-mer based approach.
 *
 * @status stable
 * @keywords streptococcus pneumoniae, serotype, k-mer, prediction, seroba
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent,conditional-logic
 * @citation seroba
 *
 * @note Database Required
 * Requires the SeroBA database to be set up using `seroba createDBs` before running.
 *
 * @input record(meta, r1, r2)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: Serotype prediction results with predicted serotype and confidence in TSV format
 *
 * @results supplemental
 * - `detailed_serogroup_info.txt`: Detailed k-mer hit information per serogroup
 */
nextflow.preview.types = true

process SEROBA_RUN {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Map,
        r1: Path,
        r2: Path
    )

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}.tsv"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.tsv"),
            files("supplemental/*", optional: true)
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
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    """
    seroba \\
        runSerotyping \\
        ${r1} ${r2} \\
        supplemental ${task.ext.args}

    # Avoid name collisions
    mv supplemental/pred.tsv ${prefix}.tsv

    # Cleanup

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seroba: \$(seroba version)
    END_VERSIONS
    """
}
