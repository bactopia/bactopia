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
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: Tab-delimited LisSero results with predicted serogroup and marker gene detection
 */
nextflow.preview.types = true

process LISSERO {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, assembly: Path): Record

    output:
    record(
        meta: meta,
        // Named field (upstream consumers access this)
        tsv: file("${prefix}.tsv"),
        // Generic fields (same convention across every module)
        results: [file("${prefix}.tsv")],
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
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${assembly} > ${assembly_name}
    fi

    lissero \\
        ${task.ext.args} \\
        ${assembly_name} \\
        > ${prefix}.tsv
    sed -i 's/^.*${assembly_name}/${assembly_name}/' ${prefix}.tsv

    # Cleanup
    rm -rf ${assembly_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lissero: \$( echo \$(lissero --version 2>&1) | sed 's/^.*LisSero //' )
    END_VERSIONS
    """
}
