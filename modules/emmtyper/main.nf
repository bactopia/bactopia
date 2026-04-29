/**
 * *emm*-typing of *Streptococcus pyogenes* (Group A Strep) assemblies.
 *
 * Uses [emmtyper](https://github.com/MDU-PHL/emmtyper) to assign *emm* types to
 * *S. pyogenes* genomes by blasting the assembly against a database of specific
 * M protein gene (*emm*) subtypes.
 *
 * @status stable
 * @keywords bacteria, streptococcus pyogenes, gas, typing, emm, virulence, m protein
 * @tags complexity:moderate input-type:single output-type:single features:database-dependent,conditional-logic
 * @citation emmtyper
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Record containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @input blastdb?
 * Path to a custom *emm* cluster BLAST database
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: Tab-delimited summary of the assigned emm type and cluster
 */
nextflow.enable.types = true

process EMMTYPER {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Record,
        fna: Path
    )
    blastdb: Path?

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
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fna} > ${fna_name}
    fi

    # Conditionally add the database if it is provided by user
    if [ "${blastdb}" == "null" ]; then
        emmtyper \\
            ${task.ext.args} \\
            ${fna_name} \\
            > ${prefix}.tsv
    else
        # Make the blast database
        makeblastdb -in ${blastdb} -dbtype nucl

        emmtyper \\
            --blast_db ${blastdb} \\
            ${task.ext.args} \\
            ${fna_name} \\
            > ${prefix}.tsv

        # Remove the blast database
        rm ${blastdb}.*
    fi

    # If 'tmp' is not in ${fna_name}, remove '.tmp' from the output files contents
    if [ ${fna_name} != *.tmp* ]; then
        sed -i 's/.tmp\t/\t/g' ${prefix}.tsv
    fi

    # Cleanup
    if [ "${is_compressed}" == "true" ]; then
        rm -rf ${fna_name}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emmtyper: \$( echo \$(emmtyper --version 2>&1) | sed 's/^.*emmtyper v//' )
    END_VERSIONS
    """
}
