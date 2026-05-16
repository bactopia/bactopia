/**
 * Filter host reads from sequencing data using minimizer-based comparison.
 *
 * Uses [deacon](https://github.com/bede/deacon) to identify and remove host reads from
 * FASTQ files using SIMD-accelerated minimizer comparison against a pre-built or custom
 * reference index. Supports paired-end, single-end, and long reads.
 *
 * @status stable
 * @keywords host, contamination, decontamination, depletion, filtering, minimizer, reads, deacon
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent,conditional-logic
 * @citation deacon
 *
 * @note Database Required
 * Requires a deacon minimizer index. Use the deacon/fetch module to download a pre-built
 * index (e.g., panhuman-1) or deacon/index to build one from a reference FASTA.
 *
 * @input record(meta, r1?, r2?, se?, lr?)
 * - `meta`: Groovy Record containing sample information
 * - `r1?`: Illumina R1 reads (paired-end forward)
 * - `r2?`: Illumina R2 reads (paired-end reverse)
 * - `se?`: Single-end Illumina reads
 * - `lr?`: Long reads (ONT/PacBio)
 *
 * @input db
 * Deacon minimizer index file (.idx) for host read filtering
 *
 * @output record(meta, special_meta, r1?, r2?, se?, lr?, scrub_report, results, logs, nf_logs, versions)
 * - `special_meta`: A simplified metadata record for downstream report joining
 * - `r1?`: Filtered paired-end forward reads
 * - `r2?`: Filtered paired-end reverse reads
 * - `se?`: Filtered single-end reads
 * - `lr?`: Filtered long reads
 * - `scrub_report`: Summary report of reads removed during filtering
 */
nextflow.enable.types = true

process DEACON_FILTER {
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
        scrub_report: file("${prefix}.scrub.report.tsv"),
        json_summary: file("${prefix}.deacon.json"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}*.scrubbed.fastq.gz"),
            files("${prefix}.scrub.report.tsv"),
            files("${prefix}.deacon.json")
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

    // Simplified meta for downstream report joining
    special_meta = record(
        name: prefix
    )

    // Pick the single-file input (se or lr)
    def single_reads = has_se ? "${se}" : "${lr}"
    if (meta.single_end) {
        """
        deacon \\
            filter \\
            --threads ${task.cpus} \\
            --summary ${prefix}.deacon.json \\
            ${task.ext.args} \\
            ${db} \\
            ${single_reads} \\
            -o ${prefix}.scrubbed.fastq.gz

        # Quick stats on reads
        zcat ${single_reads} | fastq-scan > original.json
        zcat *.scrubbed.fastq.gz | fastq-scan > scrubbed.json
        bactopia-scrubber-summary ${prefix} original.json scrubbed.json > ${prefix}.scrub.report.tsv

        # Cleanup
        rm original.json scrubbed.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            deacon: \$( deacon --version | head -n1 | sed 's/deacon //' )
            fastq-scan: \$(echo \$(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
        END_VERSIONS
        """
    } else {
        """
        deacon \\
            filter \\
            --threads ${task.cpus} \\
            --summary ${prefix}.deacon.json \\
            ${task.ext.args} \\
            ${db} \\
            ${r1} ${r2} \\
            -o ${prefix}_R1.scrubbed.fastq.gz \\
            -O ${prefix}_R2.scrubbed.fastq.gz

        # Quick stats on reads
        zcat ${r1} ${r2} | fastq-scan > original.json
        zcat *.scrubbed.fastq.gz | fastq-scan > scrubbed.json
        bactopia-scrubber-summary ${prefix} original.json scrubbed.json > ${prefix}.scrub.report.tsv

        # Cleanup
        rm original.json scrubbed.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            deacon: \$( deacon --version | head -n1 | sed 's/deacon //' )
            fastq-scan: \$(echo \$(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
        END_VERSIONS
        """
    }
}
