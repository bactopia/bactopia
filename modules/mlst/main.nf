/**
 * Automatic Multi-Locus Sequence Typing (MLST) of genome assemblies.
 *
 * Uses [mlst](https://github.com/tseemann/mlst) to scan genome assemblies against traditional
 * PubMLST schemes. It automatically detects the likely species scheme, identifies the alleles
 * for the 7 housekeeping genes, and assigns a Sequence Type (ST).
 *
 * @status stable
 * @keywords bacteria, typing, mlst, sequence type, pubmlst, alleles
 * @tags complexity:simple input-type:single output-type:single features:database-dependent,conditional-logic
 * @citation mlst
 *
 * @note Database Required
 * Requires the MLST database (derived from PubMLST) to be available.
 *
 * @input record(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @input db
 * Directory or compressed tarball containing the MLST database schemes
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: A tab-delimited summary containing the Sample, Scheme, ST, and Allele IDs
 */
nextflow.preview.types = true

process MLST {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, assembly: Path): Record
    db                : Path

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

    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    """
    # Extract database
    if [ "${is_tarball}" == "true" ]; then
        mkdir database
        tar -xzf ${db} -C database
        MLST_DB=\$(find database/ -name "mlst.fa" | sed 's=blast/mlst.fa==')
    else
        MLST_DB=\$(find ${db}/ -name "mlst.fa" | sed 's=blast/mlst.fa==')
    fi

    mlst \\
        --threads ${task.cpus} \\
        --blastdb \$MLST_DB/blast/mlst.fa \\
        --datadir \$MLST_DB/pubmlst \\
        ${task.ext.args} \\
        ${assembly} \\
        > ${prefix}.tsv

    if [[ -f "\$MLST_DB/DB_VERSION" ]]; then
        DB_VERSION=\$(cat \$MLST_DB/DB_VERSION)
    else
        DB_VERSION="custom database"
    fi

    # Cleanup
    rm -rf database/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlst: \$( echo \$(mlst --version 2>&1) | sed 's/^.*mlst //' )
        database: \$DB_VERSION
    END_VERSIONS
    """
}
