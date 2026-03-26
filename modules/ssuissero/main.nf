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
 * @input record(meta, fna)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: SsuisSero results in TSV format
 */
nextflow.preview.types = true

process SSUISSERO {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, fna: Path): Record

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
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name

    def is_compressed = fna.getName().endsWith(".gz") ? true : false
    def fna_name = fna.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fna} > ${fna_name}
    fi

    SsuisSero.sh \\
        -i ${fna_name} \\
        -o ./ \\
        -s ${prefix} \\
        -x fasta \\
        -t ${task.cpus}

    # Cleanup
    mv ${prefix}_serotyping_res.tsv ${prefix}.tsv
    rm -rf ${fna_name} blast_res/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ssuissero: ${task.ext.version}
    END_VERSIONS
    """
}
