/**
 * Finding spa types in Staphylococcus aureus.
 *
 * Uses [spaTyper](https://github.com/HCGB-IGTP/spaTyper) to determine the *spa* type of
 * *Staphylococcus aureus* genomes by identifying the repeats in the polymorphic X region
 * of the protein A gene (*spa*).
 *
 * @status stable
 * @keywords staphylococcus aureus, spa typing, repeat, mrsa, typing
 * @tags complexity:moderate input-type:single output-type:single features:compression,conditional-logic
 * @citation spatyper
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @input repeats?
 * Custom repeat sequences file
 *
 * @input repeat_order?
 * Custom repeat order file
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: spa typing results in TSV format
 */
nextflow.preview.types = true

process SPATYPER {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (meta: Map, fna: Path): Record
    repeats     : Path?
    repeat_order: Path?

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}.tsv"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.tsv")
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

    def input_args = repeats && repeat_order ? "-r ${repeats} -o ${repeat_order}" : ""
    def is_compressed = fna.getName().endsWith(".gz") ? true : false
    def fna_name = fna.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fna} > ${fna_name}
    fi

    spaTyper \\
        ${task.ext.args} \\
        ${input_args} \\
        --fasta ${fna_name} \\
        --output ${prefix}.tsv

    # Cleanup
    if [ "${is_compressed}" == "true" ]; then
        rm -rf ${fna_name}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spatyper: \$( echo \$(spaTyper --version 2>&1) | sed 's/^.*spaTyper //' )
    END_VERSIONS
    """
}
