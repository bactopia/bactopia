/**
 * In silico typing and characterization of *Bacillus cereus* group genomes.
 *
 * Uses [BTyper3](https://github.com/lmc297/BTyper3) to classify *B. cereus* group isolates.
 * It determines the PanC clade, Multi-Locus Sequence Type (MLST), and screens for virulence
 * factors, crystal toxins (Bt), and antimicrobial resistance genes.
 *
 * @status stable
 * @keywords bacteria, bacillus, cereus, typing, virulence, toxin, amr, mlst
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic
 * @citation btyper3
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: Tab-delimited Bacillus cereus group typing results including PanC clade and virulence factors
 */
nextflow.preview.types = true

process BTYPER3 {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, assembly: Path): Record

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}.tsv"),
        // Generic fields (used for publishing)
        results: files("${prefix}.tsv") + files("supplemental/*"),
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
    def assembly_name = assembly.getName().replace(".gz", "")
    """
    # Btyper3 does not accept compressed files
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${assembly} > ${assembly_name}
    fi

    btyper3 \\
        ${task.ext.args} \\
        --output ./ \\
        --input ${assembly_name}

    mv btyper3_final_results/ supplemental/
    mv supplemental/${prefix}_final_results.txt ./${prefix}.tsv

    # Cleanup
    rm -rf ${assembly_name} ${assembly_name}.njs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        btyper3: \$(echo \$(btyper3 --version 2>&1) | sed 's/^.*btyper3 //;' ))
    END_VERSIONS
    """
}
