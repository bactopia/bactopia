/**
 * Detect anti-phage defense systems using HMM profiles.
 *
 * Uses [DefenseFinder](https://github.com/mdmparis/defense-finder) to systematically search
 * protein sequences for known antiviral defense systems (e.g., CRISPR-Cas, Restriction-Modification,
 * TA systems, CBASS) using MacSyFinder and a dedicated database of HMM models.
 *
 * @status stable
 * @keywords bacteria, defense systems, antiviral, phage, crispr, restriction-modification, hmm, macsyfinder
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent,compression
 * @citation defensefinder
 *
 * @note Database Required
 * Requires the DefenseFinder HMM database to be available.
 *
 * @input record(meta, faa)
 * - `meta`: Groovy Record containing sample information
 * - `faa`: Protein sequences in FASTA format (amino acids)
 *
 * @input db
 * Directory containing the DefenseFinder models database
 *
 * @output record(meta, genes_tsv, hmmer_tsv, systems_tsv, proteins?, proteins_index?, macsydata_raw?, results, logs, nf_logs, versions)
 * - `genes_tsv`: Tab-delimited list of detected defense genes
 * - `hmmer_tsv`: Tab-delimited list of HMMER hits used for detection
 * - `systems_tsv`: Tab-delimited summary of detected defense systems
 * - `proteins?`: Protein sequences of the detected defense genes
 * - `proteins_index?`: Index file for the protein sequences
 * - `macsydata_raw?`: Compressed tarball of raw MacSyFinder data
 */
nextflow.enable.types = true

process DEFENSEFINDER_RUN {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Record,
        faa: Path
    )
    db: Path

    stage:
    stageAs faa, 'staging/faa/*'

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        genes_tsv: file("${prefix}_defense_finder_genes.tsv"),
        hmmer_tsv: file("${prefix}_defense_finder_hmmer.tsv"),
        systems_tsv: file("${prefix}_defense_finder_systems.tsv"),
        proteins: file("${prefix}.prt", optional: true),
        proteins_index: file("${prefix}.prt.idx", optional: true),
        macsydata_raw: file("${prefix}.macsydata.tar.gz", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}_defense_finder_genes.tsv"),
            files("${prefix}_defense_finder_hmmer.tsv"),
            files("${prefix}_defense_finder_systems.tsv"),
            files("${prefix}.prt", optional: true),
            files("${prefix}.prt.idx", optional: true),
            files("${prefix}.macsydata.tar.gz", optional: true)
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

    def which_cat = faa.getName().endsWith(".gz") ? "zcat" : "cat"
    """
    # Extract database
    # Use custom TMPDIR to prevent FileExistsError related to writing to same tmpdir (/tmp/tmp-macsy-cache/)
    tar -xf ${db}
    mkdir -p df-tmp/df
    TMPDIR=df-tmp/df HOME=df-tmp/ macsydata \\
        install \\
        --target defense-finder/ \\
        models/defense-finder-models-v${task.ext.defensefinder_models_version}.tar.gz

    mkdir -p df-tmp/cf
    TMPDIR=df-tmp/cf HOME=df-tmp/ macsydata \\
        install \\
        --target defense-finder/ \\
        models/CasFinder-${task.ext.casfinder_version}.tar.gz

    # DefenseFinder will attempt to gunzip the original symlink input file, so we need to
    # create a temporary uncompressed copy if the input is gzipped
    ${which_cat} ${faa} > ${prefix}.fna

    TMPDIR=df-tmp/ HOME=df-tmp/ defense-finder \\
        run \\
        ${task.ext.args} \\
        --workers ${task.cpus} \\
        --models-dir defense-finder/ \\
        ${prefix}.fna

    if [ "${task.ext.defensefinder_preserveraw}" == "true" ]; then
        tar -czf ${prefix}.macsydata.tar.gz defense-finder-tmp/
        rm -rf defense-finder-tmp/
    fi

    # Cleanup
    rm -rf models/ defense-finder/ df-tmp/ ${prefix}.fna

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        defense-finder: ${task.ext.defensefinder_version}
        defense-finder-models: ${task.ext.defensefinder_models_version}
        casfinder-models: ${task.ext.casfinder_version}
    END_VERSIONS
    """
}
