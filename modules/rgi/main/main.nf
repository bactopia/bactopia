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
 * @citation rgi
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, json, results, logs, nf_logs, versions)
 * - `tsv`: RGI results in tab-separated format
 * - `json`: RGI results in JSON format (optional)
 */
nextflow.preview.types = true

process RGI_MAIN {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, assembly: Path): Record

    output:
    record(
        meta:     meta,
        tsv:      file("${prefix}.tsv"),
        json:     file("${prefix}.json", optional: true),
        // Generic fields (used for publishing)
        results:  [file("${prefix}.tsv"), file("${prefix}.json", optional: true)],
        logs:     files("*.{log,err}", optional: true),
        nf_logs:  files(".command.*"),
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
    """
    rgi \\
        main \\
        ${task.ext.args} \\
        --clean \\
        --data wgs \\
        --num_threads ${task.cpus} \\
        --output_file ${prefix} \\
        --input_sequence ${assembly}

    # Remove empty json files
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
