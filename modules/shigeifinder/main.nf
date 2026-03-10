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
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: ShigEiFinder results in TSV format
 */
nextflow.preview.types = true

process SHIGEIFINDER {
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
        // Generic fields (used for publishing)
        results:  [file("${prefix}.tsv")],
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

    def VERSION = '1.3.2'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def is_compressed = assembly.getName().endsWith(".gz") ? true : false
    def assembly_name = assembly.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${assembly} > ${assembly_name}
    fi

    shigeifinder \\
        ${task.ext.args} \\
        --output ${prefix}.tsv \\
        -t ${task.cpus} \\
        -i ${assembly_name}

    # Cleanup
    rm -rf ${assembly_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shigeifinder: ${VERSION}
    END_VERSIONS
    """
}
