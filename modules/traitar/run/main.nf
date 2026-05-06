/**
 * Predict phenotypic traits from microbial genomes.
 *
 * Uses [Traitar](https://github.com/nick-youngblut/traitar3/) to predict phenotypic
 * traits from nucleotide sequences. Traitar annotates protein families using Pfam and
 * applies machine learning models to predict 67 diverse microbial traits.
 *
 * @status stable
 * @keywords phenotype, traits, pfam
 * @tags complexity:simple input-type:single output-type:single features:database-dependent,conditional-logic
 * @citation traitar
 *
 * @note Database Required
 * Requires a Pfam database directory (downloaded via `traitar pfam` or the download module).
 *
 * @input record(meta, fna)
 * - `meta`: Groovy Record containing sample information
 * - `fna`: Assembled contigs in FASTA format
 *
 * @input db
 * Pfam-A HMM file for Traitar
 *
 * @output record(meta, majority_tsv, single_tsv, results, logs, nf_logs, versions)
 * - `majority_tsv`: Majority-vote combined phenotype trait predictions
 * - `single_tsv`: Single-votes combined phenotype trait predictions
 */
nextflow.enable.types = true

process TRAITAR_RUN {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Record,
        fna: Path
    )
    db: Path

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        majority_tsv: file("${prefix}.majority.tsv"),
        single_tsv: file("${prefix}.single_votes.tsv"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.majority.tsv"),
            files("${prefix}.single_votes.tsv"),
            files("supplemental/*")
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    def _meta = meta
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta record
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
    # Decompress input if needed
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fna} > ${fna_name}
    else
        cp -L ${fna} ${fna_name}
    fi

    # Create input directory and sample file for traitar
    mkdir -p input_dir
    mv ${fna_name} input_dir/

    cat > samples.tsv <<-SAMPLE_EOF
    sample_file_name\tsample_name
    ${fna_name}\t${prefix}
    SAMPLE_EOF

    # Run traitar phenotype prediction
    traitar phenotype \\
        ${db} \\
        input_dir \\
        samples.tsv \\
        from_nucleotides \\
        ${prefix} \\
        -c ${task.cpus} \\
        --overwrite \\
        ${task.ext.args}

    # Rename primary output for consistency
    mkdir supplemental
    mv ${prefix}/* supplemental/
    mv supplemental/phenotype_prediction/predictions_majority-vote_combined.txt ${prefix}.majority.tsv
    mv supplemental/phenotype_prediction/predictions_single-votes_combined.txt ${prefix}.single_votes.tsv

    # Cleanup
    rm -rf input_dir/ samples.tsv ${prefix}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        traitar: \$( traitar --version 2>&1 | tail -1 )
    END_VERSIONS
    """
}
