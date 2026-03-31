/**
 * Taxonomic classification and host filtering of sequence reads.
 *
 * Uses [Kraken2](https://github.com/DerrickWood/kraken2) to assign taxonomic labels to short
 * DNA reads by examining exact k-mer matches against a large reference database. It uses the
 * Lowest Common Ancestor (LCA) algorithm to provide high-precision classification, making it
 * ideal for metagenomics or removing host contamination (scrubbing).
 *
 * Uses explicit positional record fields for reads:
 * - Input: record(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords metagenomics, taxonomy, classification, contamination, scrubbing, k-mer, lca
 * @tags complexity:complex input-type:multiple output-type:multiple features:database-dependent,conditional-logic
 * @citation kraken2
 *
 * @note Database Required
 * Requires a standard Kraken2 database (directory or tarball). Memory usage depends on database size (Standard ~50GB).
 *
 * @input record(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio) - not typically used by Kraken2
 *
 * @input db
 * Kraken2 database (Directory or compressed tarball)
 *
 * @output record(meta, special_meta, kraken2_report, scrub_report, classified, unclassified, results, logs, nf_logs, versions)
 * - `kraken2_report`: Standard Kraken2 report containing taxonomic abundance counts
 * - `scrub_report`: Summary report of reads removed during host scrubbing (optional)
 * - `special_meta`: A simplified metadata map for internal use
 * - `classified`: Reads assigned to a taxon in the database (FASTQ)
 * - `unclassified`: Reads NOT assigned to any taxon (FASTQ)
 */
nextflow.preview.types = true

process KRAKEN2 {
    tag "${prefix}"
    label 'process_high'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (meta: Map, r1: Path?, r2: Path?, se: Path?, lr: Path?): Record
    db: Path

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        special_meta: special_meta,
        kraken2_report: file("${prefix}.kraken2.report.txt"),
        scrub_report: file("${prefix}.scrub.report.tsv", optional: true),
        classified: files("*.${classified_naming}*.fastq.gz", optional: true),
        unclassified: files("*.${unclassified_naming}*.fastq.gz", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.kraken2.report.txt"),
            files("${prefix}.scrub.report.tsv", optional: true),
            files("*.${classified_naming}*.fastq.gz", optional: true),
            files("*.${unclassified_naming}*.fastq.gz", optional: true)
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    def _meta = meta
    prefix = task.ext.prefix ?: "${_meta.name}"
    output_folder = task.ext.wf == "scrubber" || task.ext.wf == "teton" ? "scrubber" : "${task.ext.process_name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope

    if (task.ext.wf == "teton") {
        meta.output_dir = "${prefix}/teton/tools/${output_folder}"
        meta.logs_dir = "${prefix}/teton/tools/${output_folder}/logs/${task.ext.logs_subdir}"
    }
    else {
        meta.output_dir = "${prefix}/tools/${output_folder}"
        meta.logs_dir = "${prefix}/tools/${output_folder}/logs/${task.ext.logs_subdir}"
    }
    meta.process_name = task.ext.process_name

    // Determine read type from explicit slots
    has_r1 = r1 != null
    has_r2 = r2 != null
    has_se = se != null
    meta.single_end = has_se && !has_r1 && !has_r2

    // Build read inputs for kraken2
    read_inputs = meta.single_end ? "${se}" : "${r1} ${r2}"

    special_meta = [:]
    special_meta.name = prefix
    def paired = meta.single_end ? "" : "--paired"
    classified_naming = task.ext.wf != "kraken2" ? "host" : "classified"
    classified = meta.single_end ? "${prefix}.${classified_naming}.fastq" : "${prefix}.${classified_naming}#.fastq"
    unclassified_naming = task.ext.wf != "kraken2" ? "scrubbed" : "unclassified"
    unclassified = meta.single_end ? "${prefix}.${unclassified_naming}.fastq" : "${prefix}.${unclassified_naming}#.fastq"
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    """
    if [ "${is_tarball}" == "true" ]; then
        mkdir database
        tar -xzf ${db} -C database
        KRAKEN_DB=\$(find database/ -name "hash.k2d" | head -1 | sed 's=hash.k2d==')
    else
        KRAKEN_DB=\$(find ${db}/ -name "hash.k2d" | head -1 | sed 's=hash.k2d==')
    fi

    k2 classify \\
        --db \$KRAKEN_DB \\
        --threads ${task.cpus} \\
        --unclassified-out ${unclassified} \\
        --classified-out ${classified} \\
        --report ${prefix}.kraken2.report.txt \\
        --output /dev/null \\
        ${paired} \\
        ${task.ext.args} \\
        ${read_inputs}

    # If scrubbing, rename and summarize
    if [ "${unclassified_naming}" == "scrubbed" ]; then
        # Rename scrubbed reads
        if [ "${meta.single_end}" == "false" ]; then
            mv ${prefix}.${unclassified_naming}_1.fastq ${prefix}_R1.scrubbed.fastq
            mv ${prefix}.${unclassified_naming}_2.fastq ${prefix}_R2.scrubbed.fastq
        fi

        # Quick stats on reads
        zcat ${read_inputs} | fastq-scan > original.json
        cat *.scrubbed.fastq | fastq-scan > scrubbed.json
        scrubber-summary.py ${prefix} original.json scrubbed.json > ${prefix}.scrub.report.tsv

        # Remove host reads and temp json files
        rm ${prefix}.host*.fastq original.json scrubbed.json
    fi

    # Cleanup database and large files produced by Kraken2
    if [ "${is_tarball}" == "true" ]; then
        rm -rf database/
    fi

    if [[ "${task.ext.keep_filtered_reads}" == "true" || "${task.ext.wf}" == "scrubber" || "${task.ext.wf}" == "teton" ]]; then
        # Compress Kraken FASTQs
        pigz -p ${task.cpus} *.fastq
    else
        # Remove filtered FASTQs
        rm *.fastq
    fi

    # Cleanup

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastq-scan: \$(echo \$(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
