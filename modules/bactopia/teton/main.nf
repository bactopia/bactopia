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
 * @input tuple(meta, classification)
 * - `meta`: Groovy Map containing sample information
 * - `classification`: Bracken species abundance report
 *
 * @output bacteria_tsv     A tab-delimited samplesheet compatible with Bactopia (--samples) for samples identified as Bacteria
 * @output nonbacteria_tsv  A tab-delimited samplesheet for samples NOT identified as Bacteria
 * @output sizemeup         A text file containing the predicted species and genome size
 * @output logs             Optional software execution logs containing warnings/errors
 * @output nf_logs          Nextflow execution scripts and logs for debugging
 * @output versions         A YAML formatted file with software versions
 */
nextflow.preview.types = true

process BACTOPIA_SAMPLESHEET {
    tag "${prefix}"
    label 'process_single'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, classification) : Tuple<Map, Set<Path>>

    output:
    bacteria_tsv    = tuple(meta, files("${prefix}.bacteria.tsv"))
    nonbacteria_tsv = tuple(meta, files("${prefix}.nonbacteria.tsv"))
    sizemeup        = tuple(meta, files("${prefix}-sizemeup.txt"))
    logs            = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs         = tuple(meta, files(".command.*"))
    versions        = tuple(meta, files("versions.yml"))

    script:
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
    teton-prepare.py \\
        ${prefix} \\
        ${prefix}-sizemeup.txt \\
        ${meta.runtype} \\
        ${meta.teton_reads} \\
        ${task.ext.outdir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sizemeup: \$(echo \$(sizemeup --version 2>&1) | sed 's/.*sizemeup-main, version //;s/ .*\$//' )
    END_VERSIONS
    """
}
