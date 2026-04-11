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
 * @input record(meta, fna)
 * - `meta`: Groovy Record containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: Tab-delimited Bacillus cereus group typing results including PanC clade and virulence factors
 *
 * @results supplemental
 * - `*_bt.txt`: Bt toxin gene screening results
 * - `*_cereulide.txt`: Cereulide synthetase gene detection results
 * - `*_amr.txt`: Antimicrobial resistance gene screening results
 * - `*_virulence.txt`: Virulence factor screening results
 * - `*_panC.txt`: PanC clade assignment details
 * - `*_mlst.txt`: Multi-locus sequence typing results
 */
nextflow.preview.types = true

process BTYPER3 {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Record,
        fna: Path
    )

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}.tsv"),
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
    def _meta = meta
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = record(
        id: "${prefix}-${task.process}",
        name: prefix,
        scope: task.ext.scope,
        output_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}",
        logs_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}",
        process_name: task.ext.process_name
    )

    def is_compressed = fna.getName().endsWith(".gz") ? true : false
    def fna_name = fna.getName().replace(".gz", "")
    """
    # Btyper3 does not accept compressed files
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fna} > ${fna_name}
    fi

    btyper3 \\
        ${task.ext.args} \\
        --output ./ \\
        --input ${fna_name}

    mv btyper3_final_results/ supplemental/
    mv supplemental/${prefix}_final_results.txt ./${prefix}.tsv

    # Cleanup
    if [ "${is_compressed}" == "true" ]; then
        rm -rf ${fna_name}
    fi
    rm -rf ${fna_name}.njs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        btyper3: \$(echo \$(btyper3 --version 2>&1) | sed 's/^.*btyper3 //;' ))
    END_VERSIONS
    """
}
