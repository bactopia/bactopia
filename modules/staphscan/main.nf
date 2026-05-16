/**
 * Genome-based surveillance analysis of Staphylococcus aureus.
 *
 * Uses [StaphSCAN](https://github.com/riccabolla/StaphSCAN) to perform genome-based
 * surveillance of *Staphylococcus aureus*, integrating species identification, MLST,
 * *spa* typing, SCCmec typing, capsular typing, and detection of virulence, biofilm,
 * and antimicrobial resistance genes.
 *
 * @status stable
 * @keywords staphylococcus aureus, surveillance, mlst, spa typing, sccmec, amr, virulence
 * @tags complexity:simple input-type:single output-type:single features:compression,conditional-logic
 * @citation staphscan
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Record containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @input db?
 * Custom MLST database directory
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: Per-sample surveillance summary with MLST, spa type, SCCmec, capsule, AGR, resistance, biofilm, and virulence results
 */
nextflow.enable.types = true

process STAPHSCAN {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Record,
        fna: Path
    )
    db: Path?

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

    // Create a new meta record
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
    def db_args = db ? "--db_mlst ${db}" : ""
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fna} > ${fna_name}
    fi

    staphscan \\
        -i ${fna_name} \\
        -o . \\
        --report ${prefix} \\
        ${db_args} ${task.ext.args}

    # Cleanup
    if [ "${is_compressed}" == "true" ]; then
        rm -rf ${fna_name}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        staphscan: \$( staphscan --version 2>&1 | sed 's/^.*staphscan //' )
    END_VERSIONS
    """
}
