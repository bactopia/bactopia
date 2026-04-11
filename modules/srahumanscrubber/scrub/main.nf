/**
 * Scrub human reads from FASTQ files.
 *
 * Uses [SRA Human Scrubber](https://github.com/ncbi/sra-human-scrubber) to identify and remove
 * human reads from sequencing data. It relies on a specific k-mer database to mask or remove
 * sequences that align to human references.
 *
 * Uses explicit positional named parameters for reads:
 * - Input: record(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords human, contamination, scrubber, decontamination, ncbi, sra
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic
 * @citation srahumanscrubber
 *
 * @input record(meta, r1?, r2?, se?, lr?)
 * - `meta`: Groovy Record containing sample information
 * - `r1?`: Illumina R1 reads (paired-end)
 * - `r2?`: Illumina R2 reads (paired-end)
 * - `se?`: Single-end Illumina reads
 * - `lr?`: Long reads (ONT/PacBio)
 *
 * @input db
 * SRA Human Scrubber database directory
 *
 * @output record(meta, special_meta, r1?, r2?, se?, lr?, scrub_report?, results, logs, nf_logs, versions)
 * - `special_meta`: A simplified metadata map for downstream report joining
 * - `r1?`: Scrubbed paired-end forward reads
 * - `r2?`: Scrubbed paired-end reverse reads
 * - `se?`: Scrubbed single-end reads
 * - `lr?`: Scrubbed long reads
 * - `scrub_report?`: Report of scrubbing statistics
 */
nextflow.preview.types = true

process SRAHUMANSCRUBBER_SCRUB {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    record (
        meta: Record,
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
        special_meta: special_meta,
        r1: file("${prefix}_R1.scrubbed.fastq.gz", optional: true),
        r2: file("${prefix}_R2.scrubbed.fastq.gz", optional: true),
        se: file("${prefix}.scrubbed.fastq.gz", optional: true),
        lr: file("${prefix}.scrubbed.fastq.gz", optional: true),
        scrub_report: file("${prefix}.scrub.report.tsv", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}*.scrubbed.fastq.gz"),
            files("${prefix}.scrub.report.tsv", optional: true)
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    def _meta = meta
    prefix = task.ext.prefix ?: "${_meta.name}"
    output_folder = task.ext.wf == "scrubber" || task.ext.wf == "teton" ? "scrubber" : "${task.ext.process_name}"

    // Determine read type from explicit slots
    has_r1 = r1 != null
    has_r2 = r2 != null
    has_se = se != null
    has_lr = lr != null

    // Create a new meta variable
    meta = record(
        id: "${prefix}-${task.process}",
        name: prefix,
        scope: task.ext.scope,
        output_dir: "${prefix}/tools/${output_folder}",
        logs_dir: "${prefix}/tools/${output_folder}/logs/${task.ext.logs_subdir}",
        process_name: task.ext.process_name,
        single_end: (has_se || has_lr) && !has_r1 && !has_r2,
        runtype: _meta.runtype != null ? _meta.runtype : (has_r1 && has_r2 ? "paired-end" : (has_lr ? "lr" : "se"))
    )
    special_meta = record(
        name: prefix
    )

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def VERSION = '2.2.1'
    if (meta.single_end) {
        """
        # Scrub human reads
        zcat ${se} | \
            scrub.sh -d ${db} -p ${task.cpus} | \
            gzip > ${prefix}.scrubbed.fastq.gz

        # Quick stats on reads
        zcat ${se} | fastq-scan > original.json
        zcat *.scrubbed.fastq.gz | fastq-scan > scrubbed.json
        bactopia-scrubber-summary ${prefix} original.json scrubbed.json > ${prefix}.scrub.report.tsv

        # Cleanup
        rm original.json scrubbed.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastq-scan: \$(echo \$(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
            sra-human-scrubber: ${VERSION}
        END_VERSIONS
        """
    } else {
        """
        # Scrub human reads
        zcat ${r1} | \
            scrub.sh -d ${db} -p ${task.cpus} | \
            gzip > ${prefix}_R1.scrubbed.fastq.gz
        zcat ${r2} | \
            scrub.sh -d ${db} -p ${task.cpus} | \
            gzip > ${prefix}_R2.scrubbed.fastq.gz

        # Quick stats on reads
        zcat ${r1} ${r2} | fastq-scan > original.json
        zcat *.scrubbed.fastq.gz | fastq-scan > scrubbed.json
        bactopia-scrubber-summary ${prefix} original.json scrubbed.json > ${prefix}.scrub.report.tsv

        # Cleanup
        rm original.json scrubbed.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastq-scan: \$(echo \$(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
            sra-human-scrubber: ${VERSION}
        END_VERSIONS
        """
    }
}
