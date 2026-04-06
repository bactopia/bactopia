/**
 * Search a protein database using a translated nucleotide query.
 *
 * Uses [BLASTX](https://blast.ncbi.nlm.nih.gov/Blast.cgi) to translate nucleotide query sequences
 * (FASTA) in all six reading frames and align them against a protein BLAST database. This is useful
 * for identifying potential coding regions in unannotated DNA.
 *
 * @status stable
 * @keywords blast, blastx, alignment, translation, protein, dna, search, fasta
 * @tags complexity:moderate input-type:multiple output-type:single features:compression
 * @citation blast
 *
 * @input record(meta, blastdb)
 * - `meta`: Groovy Map containing sample information
 * - `blastdb`: A compressed tarball containing the protein BLAST database
 *
 * @input query
 * FASTA file containing nucleotide query sequences
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: Tab-delimited translated nucleotide-to-protein alignment results (BLAST outfmt 6)
 */
nextflow.preview.types = true

process BLAST_BLASTX {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Map,
        blastdb: Path
    )
    query: Path

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}.blastx.tsv"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.blastx.tsv")
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

    def which_cat = query.getName().endsWith(".gz") ? "zcat" : "cat"
    def outcols = "sample ${task.ext.outfmt}".replace(" ", "<TAB>")
    """
    tar -xzf ${blastdb}

    ${which_cat} ${query} | \\
    blastx \\
        -num_threads ${task.cpus} \\
        -mt_mode 1 \\
        -db blastdb/${prefix}.faa \\
        -query - \\
        ${task.ext.args} \\
        -out ${prefix}.txt

    # Add column names, include column for sample name
    echo "${outcols}" | sed 's/<TAB>/\t/g' > ${prefix}.blastx.tsv
    sed 's/^/${prefix}\t/' ${prefix}.txt >> ${prefix}.blastx.tsv

    # Cleanup
    rm -rf blastdb/ ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastx: \$(blastx -version 2>&1 | sed 's/^.*blastx: //; s/ .*\$//')
    END_VERSIONS
    """
}
