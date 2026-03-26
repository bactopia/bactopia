/**
 * Identification, classification, and annotation of translated gene matches.
 *
 * Uses [GAMMA](https://github.com/rastanton/GAMMA) (Gene Allele Mutation Microbial Assessment)
 * to identify and annotate coding sequences in an assembly that match a specific gene database.
 * It is particularly useful for detecting specific targets like antimicrobial resistance genes
 * or virulence factors while accounting for potential mutations.
 *
 * @status stable
 * @keywords gene finding, annotation, homology, alignment, gamma, psl
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent,conditional-logic
 * @citation gamma
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @input db
 * The reference gene database in FASTA format
 *
 * @output record(meta, gamma, psl, gff, fasta, results, logs, nf_logs, versions)
 * - `gamma`: Main GAMMA output file containing annotated gene matches
 * - `psl`: Raw alignment details in PSL format
 * - `gff`: Gene matches in GFF3 format
 * - `fasta`: Extracted nucleotide sequences of the matched genes
 */
nextflow.preview.types = true

process GAMMA {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, fna: Path): Record
    db                   : Path

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        gamma: file("${prefix}.gamma"),
        psl: file("${prefix}.psl"),
        gff: file("${prefix}.gff", optional: true),
        fasta: file("${prefix}.fasta", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.gamma"),
            files("${prefix}.psl"),
            files("${prefix}.gff", optional: true),
            files("${prefix}.fasta", optional: true)
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
    def VERSION = '2.1'
    // Version information not provided by tool on CLI
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fna} > ${fna_name}
    fi

    GAMMA.py \\
        ${task.ext.args} \\
        ${fna_name} \\
        ${db} \\
        ${prefix}

    # Cleanup
    rm -rf ${fna_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gamma: ${VERSION}
    END_VERSIONS
    """
}
