/**
 * Detect sequence variations in the *mcr-1* colistin resistance gene.
 *
 * Uses [Mcroni](https://github.com/tseemann/mcroni) to screen genome assemblies for the
 * *mcr-1* gene (Mobilized Colistin Resistance). It extracts the gene sequence and reports
 * any variations (mutations) relative to the reference, which is critical for tracking
 * resistance to colistin, a last-resort antibiotic.
 *
 * @status stable
 * @keywords bacteria, amr, resistance, colistin, mcr-1, plasmid, variation
 * @tags complexity:simple input-type:single output-type:multiple features:conditional-logic
 * @citation mcroni
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, fa, results, logs, nf_logs, versions)
 */
nextflow.preview.types = true

process MCRONI {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, assembly: Path): Record

    output:
    record(
        meta: meta,
        tsv: file("${prefix}.tsv"),
        fa: file("${prefix}.fasta", optional: true),
        results: files("${prefix}.tsv") + files("${prefix}.fasta", optional: true),
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    def VERSION = '1.0.4'
    // Version information not provided by tool on CLI
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

    mcroni \\
        --output ${prefix} \\
        --fasta ${assembly_name}

    EX_COLS=\$(head -n 1 ${prefix}_table.tsv| tr '\\t' '\\n' | wc -l)
    OBS_COLS=\$(tail -n 1 ${prefix}_table.tsv| tr '\\t' '\\n' | wc -l)
    if [ "\$EX_COLS" != "\$OBS_COLS" ]; then
        sed -i 's/NA\$/NA\\tNA/' ${prefix}_table.tsv
    fi

    # Cleanup
    mv ${prefix}_table.tsv ${prefix}.tsv
    mv ${prefix}_sequences.fa ${prefix}.fasta
    rm -rf ${assembly_name} ${assembly_name}.ndb ${assembly_name}.not ${assembly_name}.ntf ${assembly_name}.nto

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mcroni: ${VERSION}
    END_VERSIONS
    """
}
