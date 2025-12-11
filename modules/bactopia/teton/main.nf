/**
 * Create sample sheets for Teton workflow based on species classification.
 *
 * This process executes bactopia_samplesheet to perform analysis
 *
 * @status stable
 * @keywords bactopia, sample sheet, classification, teton
 * @tags complexity:moderate input-type:single output-type:multiple
 * @citation bactopia_samplesheet
 *
 * @input tuple(meta, classification)
 * - `meta`: Groovy Map containing sample information
 * - `classification`: Classification results file
 *
 * @output bacteria_tsv    Sample sheet for bacterial samples
 * @output nonbacteria_tsv Sample sheet for non-bacterial samples
 * @output sizemeup        Genome size predictions
 * @output logs            Optional tool execution logs
 * @output nf_logs         Nextflow execution logs
 * @output versions        Software version information (YAML format)
 */
nextflow.preview.types = true

process BACTOPIA_SAMPLESHEET {
    tag "${prefix}"
    label 'process_single'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, classification) : Tuple<Map, Path>

    output:
    bacteria_tsv    = tuple(meta, file("${prefix}.bacteria.tsv"))
    nonbacteria_tsv = tuple(meta, file("${prefix}.nonbacteria.tsv"))
    sizemeup        = tuple(meta, file("${prefix}-sizemeup.txt"))
    logs            = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs         = tuple(meta, files(".command.*"))
    versions        = tuple(meta, file("versions.yml"))

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
