/**
 * Predict antibiotic resistance from assemblies.
 *
 * Uses [RGI](https://github.com/arpcard/rgi) (Resistance Gene Identifier) to predict
 * resistomes from protein or nucleotide data based on homology and SNP models using
 * the Comprehensive Antibiotic Resistance Database (CARD).
 *
 * @status stable
 * @keywords resistance, antimicrobial resistance, card, rgi, amr
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic
 * @citation rgi, diamond
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Record containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, json?, results, logs, nf_logs, versions)
 * - `tsv`: RGI results in tab-separated format
 * - `json?`: RGI results in JSON format
 */
nextflow.enable.types = true

process RGI_MAIN {
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
        json: file("${prefix}.json", optional: true),
        // Generic fields (used for publishing)
        results:  [
            files("${prefix}.tsv"),
            files("${prefix}.json", optional: true)
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
    """
    rgi \\
        main \\
        ${task.ext.args} \\
        --clean \\
        --data wgs \\
        --num_threads ${task.cpus} \\
        --output_file ${prefix} \\
        --input_sequence ${fna}

    # Cleanup
    if grep "^{}\$" ${prefix}.json; then
        rm ${prefix}.json
    fi
    mv ${prefix}.txt ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rgi: \$(rgi main --version)
    END_VERSIONS
    """
}
