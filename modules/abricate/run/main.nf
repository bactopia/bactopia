/**
 * Mass screening of contigs for antimicrobial and virulence genes.
 *
 * Screens assemblies for antimicrobial resistance and virulence genes using
 * [Abricate](https://github.com/tseemann/abricate). It bundles several databases
 * including NCBI, CARD, ResFinder, PlasmidFinder, ARG-ANNOT, and VFDB.
 *
 * @status stable
 * @keywords bacteria, assembly, fasta, antimicrobial resistance, virulence, plasmid, mobile genetic elements
 * @tags complexity:simple input-type:single output-type:single features:database-dependent
 * @citation abricate
 *
 * @note Database Included
 * Abricate bundles multiple databases including NCBI, CARD, ResFinder, PlasmidFinder,
 * ARG-ANNOT, and VFDB.
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output record(meta, report, results, logs, nf_logs, versions)
 */
nextflow.preview.types = true

process ABRICATE_RUN {
    tag "${prefix}"
    label 'process_single'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, assembly: Path): Record

    output:
    record(
        meta: meta,
        report: file("${prefix}.tsv"),
        results: [file("${prefix}.tsv")],
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
    """
    abricate \\
        ${assembly} \\
        $task.ext.args \\
        --threads $task.cpus > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' )
    END_VERSIONS
    """
}
