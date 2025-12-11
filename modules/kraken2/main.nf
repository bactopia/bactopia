/**
 * Taxonomic classification of reads using Kraken2.
 *
 * This process executes kraken2 to perform analysis
 *
 * @status stable
 * @keywords kraken2, classification, taxonomy, fastq, metagenomics
 * @tags complexity:complex input-type:multiple output-type:multiple features:archive-output, compression, conditional-logic, database-dependent, path-workarounds
 * @citation kraken2
 *
 * @note Uses EMPTY_* placeholder files for optional parameters
 * @note Requires external database to be available
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: Paired-end or single-end fastq files
 *
 * @input db
 * Kraken2 database directory or tarball
 *
 * @output kraken2_report       Kraken2 report file
 * @output scrub_report         Scrubbing report (optional)
 * @output scrub_special_report Scrub Special Report
 * @output classified           Classified reads (optional)
 * @output unclassified         Unclassified reads (optional)
 * @output classified_extra     Classified reads with extra placeholder file
 * @output unclassified_extra   Unclassified reads with extra placeholder file
 * @output logs                 Optional tool execution logs
 * @output nf_logs              Nextflow execution logs
 * @output versions             Software version information (YAML format)
 */
nextflow.preview.types = true

process KRAKEN2 {
    tag "${prefix}"
    label 'process_high'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, reads) : Tuple<Map, List<Path>>
    db             : Path

    output:
    kraken2_report       = tuple(meta, files('*.kraken2.report.txt'))
    scrub_report         = tuple(meta, files('*.scrub.report.tsv', optional: true))
    scrub_special_report = tuple(special_meta, files('*.scrub.report.tsv', optional: true))
    classified           = tuple(meta, files("*.${classified_naming}*.fastq.gz", optional: true))
    unclassified         = tuple(meta, files("*.${unclassified_naming}*.fastq.gz", optional: true))
    classified_extra     = tuple(meta, files("*.${classified_naming}*.fastq.gz", optional: true), file("EMPTY_EXTRA", optional: true))
    unclassified_extra   = tuple(meta, files("*.${unclassified_naming}*.fastq.gz", optional: true), file("EMPTY_EXTRA", optional: true))
    logs                 = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs              = tuple(meta, files(".command.*"))
    versions             = tuple(meta, file("versions.yml"))

    script:
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
    meta.single_end = reads[1] == null ? true : false
    meta.is_paired = reads[1] == null ? false : true
    meta.runtype = _meta.runtype
    special_meta = [:]
    special_meta.id = prefix
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
        KRAKEN_DB=\$(find database/ -name "hash.k2d" | sed 's=hash.k2d==')
    else
        KRAKEN_DB=\$(find ${db}/ -name "hash.k2d" | sed 's=hash.k2d==')
    fi

    kraken2 \\
        --db \$KRAKEN_DB \\
        --threads ${task.cpus} \\
        --unclassified-out ${unclassified} \\
        --classified-out ${classified} \\
        --report ${prefix}.kraken2.report.txt \\
        --gzip-compressed \\
        ${paired} \\
        ${task.ext.args} \\
        ${reads} > /dev/null

    # If scrubbing, rename and summarize
    if [ "${unclassified_naming}" == "scrubbed" ]; then
        # Rename scrubbed reads
        if [ "${meta.single_end}" == "false" ]; then
            mv ${prefix}.${unclassified_naming}_1.fastq ${prefix}_R1.scrubbed.fastq
            mv ${prefix}.${unclassified_naming}_2.fastq ${prefix}_R2.scrubbed.fastq
        fi

        # Quick stats on reads
        zcat ${reads} | fastq-scan > original.json
        cat *.scrubbed.fastq | fastq-scan > scrubbed.json
        scrubber-summary.py ${prefix} original.json scrubbed.json > ${prefix}.scrub.report.tsv

        # Remove host reads and temp json files
        rm ${prefix}.host*.fastq original.json scrubbed.json
    fi

    # Clean up database and large files produced by Kraken2
    if [ "${is_tarball}" == "true" ]; then
        rm -rf database
    fi

    if [[ "${task.ext.keep_filtered_reads}" == "true" || "${task.ext.wf}" == "scrubber" || "${task.ext.wf}" == "teton" ]]; then
        # Compress Kraken FASTQs
        pigz -p ${task.cpus} *.fastq
    else
        # Remove filtered FASTQs
        rm *.fastq
    fi

    # Used for clean-yer-reads
    touch EMPTY_EXTRA

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastq-scan: \$(echo \$(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
