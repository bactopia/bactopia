/**
 * Predict genome size and route samples based on taxonomic classification.
 *
 * Uses [SizeMeUp](https://github.com/bactopia/bactopia) to parse [Bracken](https://github.com/jenniferlu717/Bracken)
 * abundance reports, estimate the genome size for the identified species, and split samples
 * into "Bacteria" (for downstream analysis with Bactopia) and "Non-Bacteria" lists.
 *
 * @status stable
 * @keywords taxonomy, genome size, routing, filtering, bacteria, sizemeup, bracken
 * @tags complexity:simple input-type:single output-type:multiple features:conditional-logic
 * @citation bactopia, bracken, sizemeup
 *
 * @input record(meta, classification)
 * - `meta`: Groovy Map containing sample information
 * - `classification`: Bracken species abundance report
 *
 * @output record(meta, bacteria_tsv, nonbacteria_tsv, sizemeup, results, logs, nf_logs, versions)
 * - `bacteria_tsv`: A tab-delimited samplesheet compatible with Bactopia (--samples) for samples identified as Bacteria
 * - `nonbacteria_tsv`: A tab-delimited samplesheet for samples NOT identified as Bacteria
 * - `sizemeup`: A text file containing the predicted species and genome size
 */
nextflow.preview.types = true

process BACTOPIA_SAMPLESHEET {
    tag "${prefix}"
    label 'process_single'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (meta: Map, classification: Path): Record

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        bacteria_tsv: file("${prefix}.bacteria.tsv"),
        nonbacteria_tsv: file("${prefix}.nonbacteria.tsv"),
        sizemeup: file("${prefix}-sizemeup.txt"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.bacteria.tsv"),
            files("${prefix}.nonbacteria.tsv"),
            files("${prefix}-sizemeup.txt")
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
    meta.output_dir = "${prefix}/teton/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/teton/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    meta.runtype = _meta.runtype
    meta.teton_reads = _meta.teton_reads
    """
    # determine genome size and create sample sheet
    sizemeup \\
        --query ${classification} \\
        --prefix ${prefix}

    # create sample sheet
    if [ ${meta.run_type} != "ci" ]; then
        bactopia-teton-prepare \\
            ${prefix} \\
            ${prefix}-sizemeup.txt \\
            ${meta.runtype} \\
            ${meta.teton_reads} \\
            ${task.ext.outdir}
    else
        # This is a CI run, outfir path is not available
        touch ${prefix}.bacteria.tsv
        touch ${prefix}.nonbacteria.tsv
    fi

    # Cleanup

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sizemeup: \$(echo \$(sizemeup --version 2>&1) | sed 's/.*sizemeup-main, version //;s/ .*\$//' )
    END_VERSIONS
    """
}
