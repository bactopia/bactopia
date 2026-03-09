/**
 * Serotype prediction of Streptococcus suis assemblies.
 *
 * Uses [SsuisSero](https://github.com/idmc-cnr/SsuisSero) to predict the serotype of
 * *Streptococcus suis* strains from genome assemblies based on the presence of specific
 * capsular genes.
 *
 * @status stable
 * @keywords streptococcus suis, serotype, typing, prediction
 * @tags complexity:simple input-type:single output-type:single features:compression,conditional-logic
 * @citation ssuissero
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `meta`: Groovy Map containing sample information and output paths
 * - `tsv`: SsuisSero results in TSV format
 * - `results`: List of result files for publishing
 * - `logs`: Optional software execution logs containing warnings/errors
 * - `nf_logs`: Nextflow execution scripts and logs for debugging
 * - `versions`: A YAML formatted file with software versions
 */
nextflow.preview.types = true

process SSUISSERO {
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

    SsuisSero.sh \\
        -i ${assembly_name} \\
        -o ./ \\
        -s ${prefix} \\
        -x fasta \\
        -t ${task.cpus}

    # Cleanup
    mv ${prefix}_serotyping_res.tsv ${prefix}.tsv
    rm -rf ${assembly_name} blast_res/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ssuissero: ${task.ext.version}
    END_VERSIONS
    """
}
