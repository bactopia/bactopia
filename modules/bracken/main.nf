/**
 * Taxonomic classification and abundance estimation.
 *
 * Uses [Kraken2](https://github.com/DerrickWood/kraken2) to classify reads against a
 * taxonomic database, followed by [Bracken](https://github.com/jenniferlu717/Bracken)
 * (Bayesian Reestimation of Abundance with KrakEN) to estimate relative abundances at
 * a specific taxonomic level. It also generates an interactive [Krona](https://github.com/marbl/Krona/wiki)
 * plot for visualization.
 *
 * Uses explicit positional record fields for reads:
 * - Input: record(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords metagenomics, classification, taxonomy, abundance, kraken2, bracken, krona
 * @tags complexity:complex input-type:multiple output-type:multiple features:database-dependent,conditional-logic
 * @citation bracken, kraken2, krona
 *
 * @note Database Required
 * Requires a compatible Kraken2/Bracken database (tarball).
 *
 * @input record(meta, r1?, r2?, se?, lr?)
 * - `meta`: Groovy Map containing sample information
 * - `r1?`: Illumina R1 reads (paired-end)
 * - `r2?`: Illumina R2 reads (paired-end)
 * - `se?`: Single-end Illumina reads
 * - `lr?`: Long reads (ONT/PacBio) - not typically used
 *
 * @input db
 * A compressed tarball containing the Kraken2/Bracken database
 *
 * @output record(meta, tsv, special_meta, classified?, unclassified?, kraken2_report, kraken2_output?, bracken_report, krona?, abundances, classification, adjusted_abundances, results, logs, nf_logs, versions)
 * - `tsv`: Tab-delimited summary of Bracken primary and secondary species abundances
 * - `special_meta`: A simplified metadata map for internal use
 * - `classified?`: Reads classified to belong to any of the taxa on the Kraken2 database
 * - `unclassified?`: Reads not classified to belong to any of the taxa on the Kraken2 database
 * - `kraken2_report`: Kraken2 report containing stats about classified and not classified reads
 * - `kraken2_output?`: Kraken2 output file containing the taxonomic classification of each read
 * - `bracken_report`: Bracken report containing stats about classified and not classified reads
 * - `krona?`: Interactive Krona HTML visualization
 * - `abundances`: Bracken abundance estimates for each taxon
 * - `classification`: Bracken per-read classification details
 * - `adjusted_abundances`: Bracken abundance estimates adjusted for unclassified reads
 */
nextflow.preview.types = true

process BRACKEN {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Map,
        r1: Path?,
        r2: Path?,
        se: Path?,
        lr: Path?
    )
    db: Path

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}.bracken.tsv"),
        special_meta: special_meta,
        classified: files('*classified*', optional: true),
        unclassified: files('*unclassified*', optional: true),
        kraken2_report: file("${prefix}.kraken2.report.txt"),
        kraken2_output: file("${prefix}.kraken2.output.txt", optional: true),
        bracken_report: file("${prefix}.bracken.report.txt"),
        krona: files("*.krona.html", optional: true),
        abundances: file("${prefix}.bracken.abundances.txt"),
        classification: file("${prefix}.bracken.classification.txt"),
        adjusted_abundances: file("${prefix}.bracken.adjusted.abundances.txt"),
        // Generic fields (used for publishing)
        results: [
            files('*classified*', optional: true),
            files('*unclassified*', optional: true),
            files("${prefix}.bracken.tsv"),
            files("${prefix}.kraken2.report.txt"),
            files("${prefix}.kraken2.output.txt", optional: true),
            files("${prefix}.bracken.report.txt"),
            files("${prefix}.bracken.abundances.txt"),
            files("${prefix}.bracken.classification.txt"),
            files("${prefix}.bracken.adjusted.abundances.txt"),
            files("*.krona.html", optional: true),
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    def _meta = meta
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    if (task.ext.wf == "teton") {
        meta.output_dir = "${prefix}/teton/tools/${task.ext.process_name}/${task.ext.subdir}"
        meta.logs_dir = "${prefix}/teton/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    } else {
        meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
        meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    }
    meta.process_name = task.ext.process_name

    // Determine read type from explicit slots
    has_r1 = r1 != null
    has_r2 = r2 != null
    has_se = se != null
    meta.single_end = has_se && !has_r1 && !has_r2

    // Build read inputs for kraken2
    read_inputs = meta.single_end ? "${se}" : "${r1} ${r2}"
    first_read = meta.single_end ? se : r1

    special_meta = [:]
    special_meta.name = prefix

    def paired = meta.single_end ? "" : "--paired"
    classified = meta.single_end ? "${prefix}.classified.fastq" : "${prefix}.classified#.fastq"
    unclassified = meta.single_end ? "${prefix}.unclassified.fastq" : "${prefix}.unclassified#.fastq"
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    def BRACKEN_VERSION = "2.7"
    def KRAKENTOOLS_VERSION = "1.2"
    meta.teton_reads = meta.single_end ? "${se}" : "${r1},${r2}"
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
        --output ${prefix}.kraken2.output.txt \\
        ${paired} \\
        ${task.ext.args} \\
        ${read_inputs}

    # Get read length
    if [ "${task.ext.bracken_read_length}" == "0" ]; then
        OBS_READ_LENGTH=\$(zcat ${first_read} | fastq-scan -q | jq -r '.qc_stats.read_median')
        echo \$OBS_READ_LENGTH
        # Pre-built Bracken databases come with 50,75,100,150,200,250,300, split the difference
        if [ "\$OBS_READ_LENGTH" -gt 275 ]; then
            READ_LENGTH="300"
        elif [ "\$OBS_READ_LENGTH" -gt 225 ]; then
            READ_LENGTH="250"
        elif [ "\$OBS_READ_LENGTH" -gt 175 ]; then
            READ_LENGTH="200"
        elif [ "\$OBS_READ_LENGTH" -gt 125 ]; then
            READ_LENGTH="150"
        elif [ "\$OBS_READ_LENGTH" -gt 85 ]; then
            READ_LENGTH="100"
        elif [ "\$OBS_READ_LENGTH" -gt 65 ]; then
            READ_LENGTH="75"
        else
            READ_LENGTH="50"
        fi
    else
        # use user defined read length
        READ_LENGTH="${task.ext.bracken_read_length}"
    fi

    bracken \\
        ${task.ext.args2} \\
        -d \$KRAKEN_DB \\
        -r \$READ_LENGTH \\
        -i ${prefix}.kraken2.report.txt \\
        -w ${prefix}.bracken.report.txt \\
        -o bracken.temp

    # Sort bracken report by 'fraction_total_reads' (column 7)
    head -n 1 bracken.temp > ${prefix}.bracken.abundances.txt
    grep -v "fraction_total_reads\$" bracken.temp | sort -k 7 -rn -t \$'\t' >> ${prefix}.bracken.abundances.txt

    # Adjust bracken to include unclassified and produce summary
    bactopia-kraken-bracken-summary \\
        ${prefix} \\
        ${prefix}.kraken2.report.txt \\
        ${prefix}.bracken.report.txt \\
        ${prefix}.bracken.abundances.txt \\
        --max_secondary_percent ${task.ext.bracken_max_secondary_percent}

    # Create a Krona report from reports
    if [ "${task.ext.skip_krona}" == "false" ]; then
        # Kraken2
        kreport2krona.py \\
            --report ${prefix}.kraken2.report.txt \\
            --output kraken2-krona.temp
        ktImportText -o ${prefix}.kraken2.krona.html kraken2-krona.temp

        # Bracken
        kreport2krona.py \\
            --report ${prefix}.bracken.report.txt \\
            --output bracken-krona.temp
        ktImportText -o ${prefix}.bracken.krona.html bracken-krona.temp
    fi

    # Cleanup large files produced by Kraken2/Bracken
    rm *.temp
    if [ "${task.ext.kraken2_keep_raw_output}" == "false" ]; then
        # Remove kraken2 STDOUT output file
        rm ${prefix}.kraken2.output.txt
    fi

    if [ "${task.ext.kraken2_keep_filtered_reads}" == "true" ]; then
        # Compress Kraken FASTQs
        pigz -p ${task.cpus} *.fastq
    else
        # Remove filtered FASTQs
        rm *.fastq
    fi

    if [ "${is_tarball}" == "true" ]; then
        rm -rf database/
    fi

    # Cleanup

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bracken: ${BRACKEN_VERSION}
        fastq-scan: \$(echo \$(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
        jq: \$(echo \$(jq --version 2>&1) | sed 's/jq-//')
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
        krakentools: ${KRAKENTOOLS_VERSION}
        krona: \$( echo \$(ktImportText 2>&1) | sed 's/^.*KronaTools //g; s/- ktImportText.*\$//g')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
        python: \$(echo \$(python --version 2>&1) | sed 's/Python //')
    END_VERSIONS
    """
}
