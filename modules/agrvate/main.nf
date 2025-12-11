/**
 * Rapid identification of Staphylococcus aureus agr locus type and agr operon variants.
 *
 * This process executes agrvate to perform analysis
 *
 * @status stable
 * @keywords fasta, virulence, Staphylococcus aureus
 * @tags complexity:moderate input-type:single output-type:multiple features:archive-output, compression, conditional-logic
 * @citation agrvate
 *
 * @input tuple(meta, fasta)
 * - `meta`: Groovy Map containing sample information
 * - `fasta`: A Staphylococcus aureus fasta file.
 *
 * @output summary      A summary of the agrvate assessement
 * @output supplemental Supplemental
 * @output logs         Optional tool execution logs
 * @output nf_logs      Nextflow execution logs
 * @output versions     Software version information (YAML format)
 */
nextflow.preview.types = true

process AGRVATE {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, fasta) : Tuple<Map, Path>

    stage:
    stageAs 'input/*', fasta

    output:
    summary      = tuple(meta, file("${prefix}-summary.tab"))
    supplemental = tuple(meta, files("supplemental/*"))
    logs         = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs      = tuple(meta, files(".command.*"))
    versions     = tuple(meta, file("versions.yml"))

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
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = "${prefix}.fna"
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ./${fasta_name}
    else
        cat ${fasta} > ./${fasta_name}
    fi

    agrvate \\
        ${task.ext.args} \\
        -i ${fasta_name}

    mv ${meta.name}-results/ supplemental/
    mv supplemental/${meta.name}-summary.tab ./

    # Cleanup
    rm -rf ${fasta_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agrvate: \$(echo \$(agrvate -v 2>&1) | sed 's/agrvate v//;')
    END_VERSIONS
    """
}
