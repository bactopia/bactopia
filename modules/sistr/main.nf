/**
 * Serovar prediction of Salmonella assemblies.
 *
 * Uses [SISTR](https://github.com/phac-nml/sistr_cmd) (Salmonella In Silico Typing Resource) to
 * predict serovars of *Salmonella* from draft genome assemblies using core genome Multi-Locus
 * Sequence Typing (cgMLST).
 *
 * @status stable
 * @keywords salmonella, serotype, cgmlst, typing, prediction
 * @tags complexity:moderate input-type:single output-type:multiple features:compression,conditional-logic
 * @citation sistr
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @output record(meta, tsv, allele_fasta, allele_json, cgmlst_csv, results, logs, nf_logs, versions)
 * - `tsv`: SISTR prediction results in TSV format
 * - `allele_fasta`: Novel alleles in FASTA format
 * - `allele_json`: Alleles in JSON format
 * - `cgmlst_csv`: cgMLST profile in CSV format
 */
nextflow.preview.types = true

process SISTR {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (meta: Map, fna: Path): Record

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}.tsv"),
        allele_fasta: file("${prefix}-allele.fasta.gz"),
        allele_json: file("${prefix}-allele.json.gz"),
        cgmlst_csv: file("${prefix}-cgmlst.csv"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.tsv"),
            files("${prefix}-allele.fasta.gz"),
            files("${prefix}-allele.json.gz"),
            files("${prefix}-cgmlst.csv")
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

    def is_compressed = fna.getName().endsWith(".gz") ? true : false
    def fna_name = fna.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fna} > ${fna_name}
    fi

    sistr \\
        --qc \\
        ${task.ext.args} \\
        --threads ${task.cpus} \\
        --alleles-output ${prefix}-allele.json \\
        --novel-alleles ${prefix}-allele.fasta \\
        --cgmlst-profiles ${prefix}-cgmlst.csv \\
        --output-prediction ${prefix} \\
        --output-format tab \\
        ${fna_name}

    mv ${prefix}.tab ${prefix}.tsv
    gzip ${prefix}-allele.json
    gzip ${prefix}-allele.fasta

    # Cleanup
    if [ "${is_compressed}" == "true" ]; then
        rm -rf ${fna_name}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sistr: \$(echo \$(sistr --version 2>&1) | sed 's/^.*sistr_cmd //; s/ .*\$//' )
    END_VERSIONS
    """
}
