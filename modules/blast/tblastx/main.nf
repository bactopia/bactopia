/**
 * Search a translated nucleotide database using a translated nucleotide query.
 *
 * Uses [TBLASTX](https://blast.ncbi.nlm.nih.gov/Blast.cgi) to align nucleotide query sequences
 * (translated in all six frames) against a nucleotide BLAST database (also translated in all
 * six frames). This is useful for identifying distant relationships between nucleotide sequences
 * that have significant divergence but conserved protein structure.
 *
 * @status stable
 * @keywords blast, tblastx, alignment, translation, dna, search, fasta
 * @tags complexity:moderate input-type:multiple output-type:single features:compression
 * @citation blast
 *
 * @input record(meta, blastdb)
 * - `meta`: Groovy Map containing sample information
 * - `blastdb`: A compressed tarball containing the nucleotide BLAST database
 *
 * @input query
 * FASTA file containing nucleotide query sequences
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: Tab-delimited translated nucleotide-to-translated nucleotide alignment results (BLAST outfmt 6)
 */
nextflow.preview.types = true

process BLAST_TBLASTX {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, blastdb: Path): Record
    query: Path

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}.tblastx.tsv"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.tblastx.tsv")
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

    // genes -> ffn, contigs -> fna
    def db_type = task.ext.use_genes ? "ffn" : "fna"
    def which_cat = query.getName().endsWith(".gz") ? "zcat" : "cat"
    def outcols = "sample ${task.ext.outfmt}".replace(" ", "<TAB>")
    """
    tar -xzf ${blastdb}
    
    ${which_cat} ${query} | \\
    tblastx \\
        -num_threads ${task.cpus} \\
        -db blastdb/${prefix}.${db_type} \\
        -query - \\
        ${task.ext.args} \\
        -out ${prefix}.txt

    # Add column names, include column for sample name
    echo "${outcols}" | sed 's/<TAB>/\t/g' > ${prefix}.tblastx.tsv
    sed 's/^/${prefix}\t/' ${prefix}.txt >> ${prefix}.tblastx.tsv

    # Cleanup
    rm -rf blastdb/ ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tblastx: \$(tblastx -version 2>&1 | sed 's/^.*tblastx: //; s/ .*\$//')
    END_VERSIONS
    """
}
