/**
 * Identify SCCmec elements in Staphylococcus aureus genomes.
 *
 * Uses [SCCmec](https://github.com/rpetit3/sccmec) to identify the Staphylococcal Cassette
 * Chromosome mec (SCCmec) element in *Staphylococcus aureus* assemblies. It predicts the type
 * based on the presence of specific *mec* and *ccr* gene complexes.
 *
 * @status stable
 * @keywords sccmec, staphylococcus aureus, mrsa, antimicrobial resistance, typing
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic
 * @citation sccmec
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, targets, target_details, regions, regions_details, results, logs, nf_logs, versions)
 * - `meta`: Groovy Map containing sample information
 * - `tsv`: Main results file with SCCmec typing
 * - `targets`: BLAST results for target sequences
 * - `target_details`: Detailed results for target matches
 * - `regions`: BLAST results for SCCmec regions
 * - `regions_details`: Detailed results for SCCmec region matches
 * - `results`: List of all result files for publishing
 * - `logs`: Optional software execution logs containing warnings/errors
 * - `nf_logs`: Nextflow execution scripts and logs for debugging
 * - `versions`: A YAML formatted file with software versions
 */
nextflow.preview.types = true

process SCCMEC {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, assembly: Path): Record

    output:
    record(
        meta: meta,
        tsv: file("${prefix}.tsv"),
        targets: file("${prefix}.targets.blastn.tsv"),
        target_details: file("${prefix}.targets.details.tsv"),
        regions: file("${prefix}.regions.blastn.tsv"),
        regions_details: file("${prefix}.regions.details.tsv"),
        results: [file("${prefix}.tsv"), file("${prefix}.targets.blastn.tsv"), file("${prefix}.targets.details.tsv"), file("${prefix}.regions.blastn.tsv"), file("${prefix}.regions.details.tsv")],
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
    sccmec \\
        ${task.ext.args} \\
        --prefix ${prefix} \\
        --input ${assembly}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sccmec: \$(echo \$(sccmec --version 2>&1) | sed 's/.*sccmec_regions, version //;s/ .*\$//' )
        camlhmp: \$(echo \$(sccmec --version 2>&1) | sed 's/.*camlhmp, version //;s/ schema.*\$//' )
    END_VERSIONS
    """
}
