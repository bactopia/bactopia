/**
 * Search a protein database using a protein query.
 *
 * Uses [BLASTP](https://blast.ncbi.nlm.nih.gov/Blast.cgi) to align amino acid query sequences
 * (FASTA) against a protein BLAST database. It is used to identify homologous proteins.
 *
 * @status stable
 * @keywords blast, blastp, alignment, protein, amino acid, search, fasta
 * @tags complexity:moderate input-type:multiple output-type:single features:compression
 * @citation blast
 *
 * @input tuple(meta, blastdb)
 * - `meta`: Groovy Map containing sample information
 * - `blastdb`: A compressed tarball containing the protein BLAST database
 *
 * @input query
 * FASTA file containing amino acid query sequences
 *
 * @output tsv       A tab-delimited summary of alignments (standard BLAST outfmt 6)
 * @output logs      Optional software execution logs containing warnings/errors
 * @output nf_logs   Nextflow execution scripts and logs for debugging
 * @output versions  A YAML formatted file with software versions
 */
nextflow.preview.types = true

process BLAST_BLASTP {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, blastdb) : Tuple<Map, Path>
    query            : Path

    output:
    tsv      = tuple(meta, files('*.blastp.tsv'))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, file("versions.yml"))

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
    def which_cat = query.getName().endsWith(".gz") ? "zcat" : "cat"
    def outcols = "sample ${task.ext.outfmt}".replace(" ", "<TAB>")
    """
    tar -xzf ${blastdb}
    
    ${which_cat} ${query} | \\
    blastp \\
        -num_threads ${task.cpus} \\
        -mt_mode 1 \\
        -db blastdb/${prefix}.faa \\
        -query - \\
        ${task.ext.args} \\
        -out ${prefix}.txt

    # Add column names, include column for sample name
    echo "${outcols}" | sed 's/<TAB>/\t/g' > ${prefix}.blastp.tsv
    sed 's/^/${prefix}\t/' ${prefix}.txt >> ${prefix}.blastp.tsv

    # Cleanup
    rm -rf blastdb/ ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastp: \$(blastp -version 2>&1 | sed 's/^.*blastp: //; s/ .*\$//')
    END_VERSIONS
    """
}
