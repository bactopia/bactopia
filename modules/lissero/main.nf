/**
 * Predict *Listeria monocytogenes* serogroup.
 *
 * Uses [LisSero](https://github.com/MDU-PHL/LisSero) to predict the serogroup of
 * *L. monocytogenes* isolates. It simulates a PCR assay by detecting specific marker genes
 * (lmo1118, lmo0737, ORF2110, ORF2819, prs) to assign the isolate to one of the major
 * molecular serogroups (IIa, IIb, IIc, IVb).
 *
 * @status stable
 * @keywords bacteria, listeria, monocytogenes, serotype, serogroup, typing, pcr
 * @tags complexity:simple input-type:single output-type:single features:compression,conditional-logic
 * @citation lissero
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Record containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: Tab-delimited LisSero results with predicted serogroup and marker gene detection
 */
nextflow.enable.types = true

process LISSERO {
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
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fna} > ${fna_name}
    fi

    lissero \\
        ${task.ext.args} \\
        ${fna_name} \\
        > ${prefix}.tsv
    sed -i 's/^.*${fna_name}/${fna_name}/' ${prefix}.tsv

    # Cleanup
    if [ "${is_compressed}" == "true" ]; then
        rm -rf ${fna_name}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lissero: \$( echo \$(lissero --version 2>&1) | sed 's/^.*LisSero //' )
    END_VERSIONS
    """
}
