/**
 * Genotyping and screening of *Klebsiella* genome assemblies.
 *
 * Uses [Kleborate](https://github.com/katholt/Kleborate) to screen *Klebsiella* assemblies
 * for Multi-Locus Sequence Type (MLST), species identity, antimicrobial resistance determinants,
 * virulence plasmids (e.g., *ybt*, *iuc*, *iro*), and capsular serotype prediction (K and O loci).
 *
 * @status stable
 * @keywords bacteria, klebsiella, amr, virulence, typing, mlst, serotype, k-locus, o-locus
 * @tags complexity:moderate input-type:single output-type:single features:database-dependent,conditional-logic
 * @citation kleborate, kaptive
 *
 * @note Database Bundled
 * Kleborate bundles the required databases for species identification, MLST,
 * and virulence/resistance gene detection.
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Record containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: Tab-delimited Kleborate results with species, MLST, virulence, and resistance predictions
 */
nextflow.preview.types = true

process KLEBORATE {
    tag "${prefix}"
    label 'process_medium'

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
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.tsv")
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
    mkdir results/
    kleborate \\
        ${task.ext.args} \\
        --outdir results/ \\
        --assemblies ${fna}

    # Rename output file to include the prefix name
    find results/ -name "*output.txt" -print0 | while read -d \$'\0' file; do mv "\$file" "${prefix}.txt"; done

    # Negative results will not have an output file
    if [ ! -f "${prefix}.txt" ]; then
        touch "${prefix}.txt"
    fi

    # Cleanup
    mv ${prefix}.txt ${prefix}.tsv
    rm -rf results/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kleborate: \$( echo \$(kleborate --version | sed 's/Kleborate v//;'))
    END_VERSIONS
    """
}
