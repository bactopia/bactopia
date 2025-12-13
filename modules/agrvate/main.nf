/**
 * Determine the agr locus type and operon variants in Staphylococcus aureus.
 *
 * Uses [AgrVATE](https://github.com/VishnuRaghuram94/AgrVATE) to type the accessory gene
 * regulator (agr) locus, a quorum sensing system critical for *Staphylococcus aureus* virulence.
 *
 * @status stable
 * @keywords bacteria, assembly, fasta, typing, virulence, staphylococcus, aureus, agr
 * @tags complexity:moderate input-type:single output-type:multiple features:compression,conditional-logic
 * @citation agrvate
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled Staphylococcus aureus contigs in FASTA format
 *
 * @output summary      A tab-delimited report containing the assigned agr type and additional details
 * @output supplemental Supplemental output files from [AgrVATE](https://github.com/VishnuRaghuram94/AgrVATE)
 * @output logs         Optional software execution logs containing warnings/errors
 * @output nf_logs      Nextflow execution scripts and logs for debugging
 * @output versions     A YAML formatted file with software versions
 */
nextflow.preview.types = true

process AGRVATE {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, assembly) : Tuple<Map, Path>

    stage:
    stageAs 'input/*', assembly

    output:
    summary      = tuple(meta, file("${prefix}.tsv"))
    supplemental = tuple(meta, files("supplemental/*"))
    logs         = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs      = tuple(meta, files(".command.*"))
    versions     = tuple(meta, files("versions.yml"))

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

    def is_compressed = assembly.getName().endsWith(".gz") ? true : false
    def assembly_name = "${prefix}.fna"
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${assembly} > ./${assembly_name}
    else
        cat ${assembly} > ./${assembly_name}
    fi

    agrvate \\
        ${task.ext.args} \\
        -i ${assembly_name}

    mv ${prefix}-results/ supplemental/
    mv supplemental/${prefix}-summary.tab ./${prefix}.tsv

    # Cleanup
    rm -rf ${assembly_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agrvate: \$(echo \$(agrvate -v 2>&1) | sed 's/agrvate v//;')
    END_VERSIONS
    """
}
