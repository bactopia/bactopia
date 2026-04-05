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
 * @input record(meta, fna)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled Staphylococcus aureus contigs in FASTA format
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: Tab-delimited summary of agr locus type and operon variants
 *
 * @results supplemental
 * - `*-agr_gp.tab`: agr group allele BLAST results
 * - `*-agr_operon.fna`: Extracted agr operon nucleotide sequence
 * - `*-frameshifts.tab`: Detected frameshift mutations in the agr operon
 */
nextflow.preview.types = true

process AGRVATE {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (meta: Map, fna: Path): Record

    stage:
    stageAs "staging/fna/*", fna

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}.tsv"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.tsv"),
            files("${prefix}.tab"),
            files("supplemental/*")
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

    def is_compressed = fna.getName().endsWith(".gz") ? true : false
    def fna_name = fna.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fna} > ./${fna_name}
    else
        # agrvate does not support symlinks
        cp ${fna} ./${fna_name}
    fi

    agrvate \\
        ${task.ext.args} \\
        -i ${fna_name}

    mv ${prefix}-results/ supplemental/
    mv *.tab supplemental/
    mv supplemental/${prefix}-summary.tab ./${prefix}.tsv

    # Cleanup
    if [ "${is_compressed}" == "true" ]; then
        rm -rf ${fna_name}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agrvate: \$(echo \$(agrvate -v 2>&1) | sed 's/agrvate v//;')
    END_VERSIONS
    """
}
