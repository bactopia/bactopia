/**
 * Predict O-antigen serogroup of Pseudomonas aeruginosa isolates.
 *
 * Uses [Pasty](https://github.com/rpetit3/pasty) (in silico serogrouping of *Pseudomonas aeruginosa* isolates)
 * to predict the O-antigen serogroup by searching the genome assembly for specific serogroup-associated
 * genes within the O-antigen locus.
 *
 * @status stable
 * @keywords bacteria, pseudomonas aeruginosa, serogroup, o-antigen, typing, blast
 * @tags complexity:simple input-type:single output-type:multiple features:conditional-logic
 * @citation pasty
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, blast, details, results, logs, nf_logs, versions)
 * - `tsv`: A tab-delimited summary file with the predicted O-antigen serogroup
 * - `blast`: A tab-delimited file of all raw BLAST hits used for the prediction
 * - `details`: A tab-delimited file with detailed gene hits for each serogroup tested
 */
nextflow.preview.types = true

process PASTY {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (meta: Map, fna: Path): Record

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}.tsv"),
        blast: file("${prefix}.blastn.tsv"),
        details: file("${prefix}.details.tsv"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.tsv"),
            files("${prefix}.blastn.tsv"),
            files("${prefix}.details.tsv")
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
    pasty \\
        ${task.ext.args} \\
        --prefix ${prefix} \\
        --input ${fna}

    # Cleanup

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pasty: \$(echo \$(pasty --version 2>&1) | sed 's/.*pasty, version //;s/ .*\$//' )
        camlhmp: \$(echo \$(pasty --version 2>&1) | sed 's/.*camlhmp, version //;s/ schema.*\$//' )
    END_VERSIONS
    """
}
