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
 * @input record(meta, fna)
 * - `meta`: Groovy Record containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, targets, target_details, regions, regions_details, results, logs, nf_logs, versions)
 * - `tsv`: Main results file with SCCmec typing
 * - `targets`: BLAST results for target sequences
 * - `target_details`: Detailed results for target matches
 * - `regions`: BLAST results for SCCmec regions
 * - `regions_details`: Detailed results for SCCmec region matches
 */
nextflow.preview.types = true

process SCCMEC {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Record,
        fna: Path
    )

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}.tsv"),
        targets: file("${prefix}.targets.blastn.tsv"),
        target_details: file("${prefix}.targets.details.tsv"),
        regions: file("${prefix}.regions.blastn.tsv"),
        regions_details: file("${prefix}.regions.details.tsv"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.tsv"),
            files("${prefix}.targets.blastn.tsv"),
            files("${prefix}.targets.details.tsv"),
            files("${prefix}.regions.blastn.tsv"),
            files("${prefix}.regions.details.tsv")
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
    """
    sccmec \\
        ${task.ext.args} \\
        --prefix ${prefix} \\
        --input ${fna}

    # Cleanup

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sccmec: \$(echo \$(sccmec --version 2>&1) | sed 's/.*sccmec_regions, version //;s/ .*\$//' )
        camlhmp: \$(echo \$(sccmec --version 2>&1) | sed 's/.*camlhmp, version //;s/ schema.*\$//' )
    END_VERSIONS
    """
}
