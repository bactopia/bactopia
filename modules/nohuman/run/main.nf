/**
 * Remove human reads from sequencing data.
 *
 * Uses [nohuman](https://github.com/mbhall88/nohuman) to classify and remove human reads
 * from FASTQ files using a Kraken2 database built from Human Pangenome Reference Consortium
 * (HPRC) genomes. Supports paired-end and single-end Illumina reads.
 *
 * @status stable
 * @keywords human, contamination, decontamination, scrubbing, reads, kraken2, nohuman
 * @tags complexity:moderate input-type:single output-type:multiple features:database-dependent,conditional-logic
 * @citation kraken2
 *
 * @note Database Required
 * Requires the nohuman Kraken2 database. Use the nohuman/download module or
 * provide a pre-existing database via --nohuman_db.
 *
 * @input record(meta, r1?, r2?, se?, lr?)
 * - `meta`: Groovy Record containing sample information
 * - `r1?`: Illumina R1 reads (paired-end forward)
 * - `r2?`: Illumina R2 reads (paired-end reverse)
 * - `se?`: Single-end Illumina reads
 * - `lr?`: Long reads (ONT/PacBio)
 *
 * @input db
 * Directory or compressed tarball containing the nohuman Kraken2 database
 *
 * @output record(meta, special_meta, r1?, r2?, se?, lr?, scrub_report, results, logs, nf_logs, versions)
 * - `special_meta`: A simplified metadata map for downstream report joining
 * - `r1?`: Scrubbed paired-end forward reads
 * - `r2?`: Scrubbed paired-end reverse reads
 * - `se?`: Scrubbed single-end reads
 * - `lr?`: Scrubbed long reads
 * - `scrub_report`: Summary report of reads removed during scrubbing
 */
nextflow.preview.types = true

process NOHUMAN_RUN {
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
        // Generic fields (used for publishing)
        results: [
            files("${prefix}*.scrubbed.fastq.gz"),
            files("${prefix}.scrub.report.tsv")
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
    has_lr = lr != null

    // Create a new meta variable
    meta = record(
        id: "${prefix}-${task.process}",
        name: prefix,
        scope: task.ext.scope,
        output_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}",
        logs_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}",
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
    def is_tarball = db.getName().endsWith(".tar.gz") ? true : false
    def output_type = "--output-type g"
    if (meta.single_end) {
        """
        if [ "${is_tarball}" == "true" ]; then
            mkdir database
            tar -xzf ${db} -C database
            DB_PATH=\$(find database/ -name "hash.k2d" | head -1 | sed 's/hash.k2d//')
        else
            DB_PATH=\$(find ${db}/ -name "hash.k2d" | head -1 | sed 's/hash.k2d//')
        fi

        nohuman \\
            --db \$DB_PATH \\
            --threads ${task.cpus} \\
            ${output_type} \\
            ${task.ext.args} \\
            --out1 ${prefix}.scrubbed.fastq.gz \\
            ${single_reads}

        # Quick stats on reads
        zcat ${single_reads} | fastq-scan > original.json
        zcat *.scrubbed.fastq.gz | fastq-scan > scrubbed.json
        bactopia-scrubber-summary ${prefix} original.json scrubbed.json > ${prefix}.scrub.report.tsv

        # Cleanup
        if [ "${is_tarball}" == "true" ]; then
            rm -rf database/
        fi
        rm original.json scrubbed.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastq-scan: \$(echo \$(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
            nohuman: \$( nohuman --version 2>&1 | sed 's/nohuman //' )
        END_VERSIONS
        """
    } else {
        """
        if [ "${is_tarball}" == "true" ]; then
            mkdir database
            tar -xzf ${db} -C database
            DB_PATH=\$(find database/ -name "hash.k2d" | head -1 | sed 's/hash.k2d//')
        else
            DB_PATH=\$(find ${db}/ -name "hash.k2d" | head -1 | sed 's/hash.k2d//')
        fi

        nohuman \\
            --db \$DB_PATH \\
            --threads ${task.cpus} \\
            ${output_type} \\
            ${task.ext.args} \\
            --out1 ${prefix}_R1.scrubbed.fastq.gz \\
            --out2 ${prefix}_R2.scrubbed.fastq.gz \\
            ${r1} ${r2}

        # Quick stats on reads
        zcat ${r1} ${r2} | fastq-scan > original.json
        zcat *.scrubbed.fastq.gz | fastq-scan > scrubbed.json
        bactopia-scrubber-summary ${prefix} original.json scrubbed.json > ${prefix}.scrub.report.tsv

        # Cleanup
        if [ "${is_tarball}" == "true" ]; then
            rm -rf database/
        fi
        rm original.json scrubbed.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastq-scan: \$(echo \$(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
            nohuman: \$( nohuman --version 2>&1 | sed 's/nohuman //' )
        END_VERSIONS
        """
    }
}
