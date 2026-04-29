/**
 * Shigella and EIEC serotyping from assemblies.
 *
 * Uses [ShigEiFinder](https://github.com/LanLab/ShigEiFinder) to differentiate *Shigella* and
 * Enteroinvasive *E. coli* (EIEC) and predict their serotypes from genome assemblies. It utilizes
 * cluster-specific marker genes to distinguish these closely related pathovars.
 *
 * @status stable
 * @keywords shigella, eiec, serotype, identification, cluster, virulence
 * @tags complexity:simple input-type:single output-type:single features:conditional-logic
 * @citation shigeifinder
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Record containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: ShigEiFinder results in TSV format
 */
nextflow.enable.types = true

process SHIGEIFINDER {
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
        results:  [
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

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def VERSION = '1.3.2'
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fna} > ${fna_name}
    fi

    shigeifinder \\
        ${task.ext.args} \\
        --output ${prefix}.tsv \\
        -t ${task.cpus} \\
        -i ${fna_name}

    # Cleanup
    if [ "${is_compressed}" == "true" ]; then
        rm -rf ${fna_name}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shigeifinder: ${VERSION}
    END_VERSIONS
    """
}
