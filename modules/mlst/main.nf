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
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @input db
 * Directory or compressed tarball containing the MLST database schemes
 *
 * @output tsv       A tab-delimited summary containing the Sample, Scheme, ST, and Allele IDs
 * @output logs      Optional software execution logs containing warnings/errors
 * @output nf_logs   Nextflow execution scripts and logs for debugging
 * @output versions  A YAML formatted file with software versions
 */
nextflow.preview.types = true

process MLST {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, assembly) : Tuple<Map, Set<Path>>
    db             : Path

    output:
    tsv      = tuple(meta, files("${prefix}.tsv"))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, files("versions.yml"))

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
