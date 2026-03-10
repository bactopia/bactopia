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
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled Staphylococcus aureus contigs in FASTA format
 *
 * @output record(meta, summary, results, logs, nf_logs, versions)
 * - `summary`: Tab-delimited summary of agr locus type and operon variants
 */
nextflow.preview.types = true

process AGRVATE {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, assembly: Path): Record

    stage:
    stageAs 'input/*', assembly

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        summary: file("${prefix}.tsv"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.tsv"),
            files("supplemental/*")
        ],
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

    def is_compressed = assembly.getName().endsWith(".gz") ? true : false
    def assembly_name = "${prefix}.fna"
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${assembly} > ./${assembly_name}
    else
        cp ${assembly} ./${assembly_name}
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
