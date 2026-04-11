/**
 * Estimate bacterial species abundance from metagenomic reads.
 *
 * Uses [MIDAS](https://github.com/snayfach/MIDAS) (Metagenomic Intra-Species Diversity Analysis System)
 * to estimate the abundance of bacterial species in metagenomic data. It maps reads to a database
 * of universal single-copy marker genes (15 genes) to provide accurate coverage and relative
 * abundance estimates.
 *
 * Uses explicit positional record fields for reads:
 * - Input: record(meta, r1, r2, se) where each read slot is Path?
 *
 * @status stable
 * @keywords metagenomics, abundance, species, midas, marker genes, diversity
 * @tags complexity:moderate input-type:multiple output-type:multiple features:database-dependent,conditional-logic
 * @citation midas
 *
 * @note Database Required
 * Requires a compatible MIDAS database (containing marker gene sequences and taxonomy).
 *
 * @input record(meta, r1?, r2?, se?)
 * - `meta`: Groovy Record containing sample information
 * - `r1?`: Illumina R1 reads (paired-end)
 * - `r2?`: Illumina R2 reads (paired-end)
 * - `se?`: Single-end Illumina reads
 *
 * @input db
 * Directory containing the MIDAS database
 *
 * @output record(meta, tsv, abundances, adjusted_abundances, results, logs, nf_logs, versions)
 * - `tsv`: A tab-delimited summary of species abundance and coverage
 * - `abundances`: Detailed species abundance profile (*.abundances.txt)
 * - `adjusted_abundances`: Relative abundance estimates adjusted for genome size (*.adjusted.abundances.txt)
 */
nextflow.preview.types = true
nextflow.enable.moduleBinaries = true

process MIDAS_SPECIES {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Record,
        r1: Path?,
        r2: Path?,
        se: Path?
    )
    db: Path

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}.midas.tsv"),
        abundances: file("${prefix}.midas.abundances.txt"),
        adjusted_abundances: file("${prefix}.midas.adjusted.abundances.txt"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.midas.tsv"),
            files("${prefix}.midas.abundances.txt"),
            files("${prefix}.midas.adjusted.abundances.txt")
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    def _meta = meta
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Determine read type from explicit slots
    has_r1 = r1 != null
    has_r2 = r2 != null
    has_se = se != null

    // Create a new meta variable
    meta = record(
        id: "${prefix}-${task.process}",
        name: prefix,
        scope: task.ext.scope,
        output_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}",
        logs_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}",
        process_name: task.ext.process_name,
        single_end: has_se && !has_r1 && !has_r2
    )

    def read_opts = meta.single_end ? "-1 ${se}" : "-1 ${r1} -2 ${r2}"
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def VERSION = '1.3.2'
    """
    if [ "${is_tarball}" == "true" ]; then
        mkdir database
        tar -xzf ${db} -C database
        MIDAS_DB=\$(find database/ -name "genome_info.txt" | sed 's=genome_info.txt==')
    else
        MIDAS_DB=\$(find ${db}/ -name "genome_info.txt" | sed 's=genome_info.txt==')
    fi

    run_midas.py \\
        species \\
        results \\
        ${read_opts} \\
        ${task.ext.args} \\
        -d \${MIDAS_DB} \\
        -t ${task.cpus}

    mv results/species/species_profile.txt ${prefix}.midas.abundances.txt
    midas-summary.py ${prefix} ${prefix}.midas.abundances.txt

    # Cleanup
    rm -rf results/
    if [ "${is_tarball}" == "true" ]; then
        rm -rf database/
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        midas: ${VERSION}
    END_VERSIONS
    """
}
