/**
 * Assess genome quality using machine learning.
 *
 * Uses [CheckM2](https://github.com/chklovski/CheckM2) to predict the completeness and
 * contamination of genome assemblies. Unlike the original CheckM, it uses a gradient boost
 * machine learning model to predict quality without relying on lineage-specific marker sets,
 * making it more accurate for novel or reduced genomes.
 *
 * @status stable
 * @keywords quality control, completeness, contamination, machine learning, bacteria, archaea
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent,conditional-logic
 * @citation checkm2
 *
 * @note Database Required
 * Requires the CheckM2 database (Diamond database file) to be available.
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Map containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @input db
 * The CheckM2 database file (*.dmnd)
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: Tab-delimited report of quality metrics (Completeness, Contamination)
 *
 * @results supplemental
 * - `protein_files/*.faa.gz`: Predicted protein sequences from Prodigal (compressed)
 * - `diamond_output/`: DIAMOND alignment results against the CheckM2 database
 */
nextflow.preview.types = true

process CHECKM2_PREDICT {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, fna: Path): Record
    db: Path

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}.tsv"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.tsv"),
            files("supplemental/*")
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
    # Decompress fasta file if compressed
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fna} > ${fna_name}
    fi

    # Check if db is a directory - if so, find the diamond database
    if [ -d "${db}" ]; then
        CHECKM2_DB=\$(find ${db}/ -name "*.dmnd")
    else
        CHECKM2_DB=${db}
    fi

    checkm2 \\
        predict \\
        --output-directory supplemental \\
        --threads ${task.cpus} \\
        --database_path \$CHECKM2_DB \\
        ${task.ext.args} \\
        --input ${fna_name}

    mv supplemental/checkm2.log ./
    mv supplemental/quality_report.tsv ./${prefix}.tsv

    # Cleanup
    gzip supplemental/protein_files/*.faa
    if [ "${is_compressed}" == "true" ]; then
        rm -rf ${fna_name}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$(checkm2 --version)
    END_VERSIONS
    """
}
