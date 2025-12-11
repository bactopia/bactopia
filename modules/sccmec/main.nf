/**
 * Identify SCCmec elements in Staphylococcus aureus genomes.
 *
 * This process executes sccmec to perform analysis
 *
 * @status stable
 * @keywords SCCmec, Staphylococcus aureus, antimicrobial resistance, cassette chromosome
 * @tags complexity:moderate input-type:single output-type:multiple
 * @citation sccmec
 *
 * @input tuple(meta, fasta)
 * - `meta`: Groovy Map containing sample information
 * - `fasta`: FASTA file containing the genome assembly
 *
 * @output tsv             Main results file with SCCmec typing
 * @output targets         BLAST results for target sequences
 * @output target_details  Detailed results for target matches
 * @output regions         BLAST results for SCCmec regions
 * @output regions_details Detailed results for SCCmec region matches
 * @output logs            Optional tool execution logs
 * @output nf_logs         Nextflow execution logs
 * @output versions        Software version information (YAML format)
 */
nextflow.preview.types = true

process SCCMEC {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta) : Tuple<Map, Path>

    output:
    tsv             = tuple(meta, file("${prefix}.tsv"))
    targets         = tuple(meta, file("${prefix}.targets.blastn.tsv"))
    target_details  = tuple(meta, file("${prefix}.targets.details.tsv"))
    regions         = tuple(meta, file("${prefix}.regions.blastn.tsv"))
    regions_details = tuple(meta, file("${prefix}.regions.details.tsv"))
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
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    """
    sccmec \\
        ${task.ext.args} \\
        --prefix ${prefix} \\
        --input ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sccmec: \$(echo \$(sccmec --version 2>&1) | sed 's/.*sccmec_regions, version //;s/ .*\$//' )
        camlhmp: \$(echo \$(sccmec --version 2>&1) | sed 's/.*camlhmp, version //;s/ schema.*\$//' )
    END_VERSIONS
    """
}
