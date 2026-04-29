/**
 * Search a translated nucleotide database using a protein query.
 *
 * Uses [TBLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi) to align amino acid query sequences
 * (FASTA) against a nucleotide BLAST database that has been dynamically translated in all six
 * reading frames. This is useful for finding gene homologs in unannotated genomic DNA.
 *
 * @status stable
 * @keywords blast, tblastn, alignment, translation, protein, dna, search, fasta
 * @tags complexity:moderate input-type:multiple output-type:single features:compression
 * @citation blast
 *
 * @input record(meta, blastdb)
 * - `meta`: Groovy Record containing sample information
 * - `blastdb`: A compressed tarball containing the nucleotide BLAST database
 *
 * @input query
 * FASTA file containing amino acid query sequences
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: Tab-delimited protein-to-translated nucleotide alignment results (BLAST outfmt 6)
 */
nextflow.enable.types = true

process BLAST_TBLASTN {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Record,
        blastdb: Path
    )
    query: Path

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}.tblastn.tsv"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.tblastn.tsv")
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

    // genes -> ffn, contigs -> fna
    def db_type = task.ext.use_genes ? "ffn" : "fna"
    def which_cat = query.getName().endsWith(".gz") ? "zcat" : "cat"
    def outcols = "sample ${task.ext.outfmt}".replace(" ", "<TAB>")
    """
    tar -xzf ${blastdb}

    ${which_cat} ${query} | \\
    tblastn \\
        -num_threads ${task.cpus} \\
        -mt_mode 1 \\
        -db blastdb/${prefix}.${db_type} \\
        -query - \\
        ${task.ext.args} \\
        -out ${prefix}.txt

    # Add column names, include column for sample name
    echo "${outcols}" | sed 's/<TAB>/\t/g' > ${prefix}.tblastn.tsv
    sed 's/^/${prefix}\t/' ${prefix}.txt >> ${prefix}.tblastn.tsv

    # Cleanup
    rm -rf blastdb/ ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tblastn: \$(tblastn -version 2>&1 | sed 's/^.*tblastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
