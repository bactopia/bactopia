/**
 * Create a SNP distance matrix from a multiple sequence alignment.
 *
 * Uses [snp-dists](https://github.com/tseemann/snp-dists) to read a FASTA alignment and
 * compute a pairwise SNP distance matrix between all sequences.
 *
 * @status stable
 * @keywords snp, distance, matrix, alignment, phylogeny
 * @tags complexity:simple input-type:single output-type:single features:conditional-logic
 * @citation snpdists
 *
 * @input record(meta, alignment)
 * - `meta`: Groovy Map containing sample information
 * - `alignment`: Multiple sequence alignment in FASTA format
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `meta`: Groovy Map containing sample information and output paths
 * - `tsv`: Pairwise SNP distance matrix in TSV format
 * - `results`: List of result files for publishing
 * - `logs`: Optional software execution logs containing warnings/errors
 * - `nf_logs`: Nextflow execution scripts and logs for debugging
 * - `versions`: A YAML formatted file with software versions
 */
nextflow.preview.types = true

process SNPDISTS {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, msa: Path): Record

    output:
    record(
        meta: meta,
        tsv: file("${prefix}.tsv"),
        results: [file("${prefix}.tsv")],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"
    process_name = _meta.process_name ?: task.ext.process_name

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.ext.process_name}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = ""
    meta.logs_dir = "${process_name}/logs/"
    meta.process_name = process_name
    """
    snp-dists \\
        ${task.ext.args} \\
        ${msa} > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snp-dists: \$(snp-dists -v 2>&1 | sed 's/snp-dists //;')
    END_VERSIONS
    """
}
