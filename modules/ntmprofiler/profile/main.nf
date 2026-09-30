/**
 * Identify NTM species and predict antimicrobial resistance from sequencing reads.
 *
 * Uses [NTM-Profiler](https://github.com/jodyphelan/NTM-Profiler) to identify
 * nontuberculous Mycobacterium species and, where a compatible database exists,
 * detect resistance-associated variants and lineage barcodes.
 *
 * @status stable
 * @keywords mycobacterium, ntm, species identification, antimicrobial resistance, variant calling
 * @tags complexity:moderate input-type:multiple output-type:multiple features:database-dependent,conditional-input
 * @citation ntmprofiler
 *
 * @note Database Selection
 * Uses a database prepared once by the NTM-Profiler update module unless a
 * custom database is provided.
 *
 * @input record(meta, r1?, r2?, se?, lr?)
 * - `meta`: Groovy Record containing sample information
 * - `r1?`: Illumina R1 reads (paired-end forward)
 * - `r2?`: Illumina R2 reads (paired-end reverse)
 * - `se?`: Single-end Illumina reads
 * - `lr?`: Long reads (ONT/PacBio)
 *
 * @input db
 * Directory or compressed tarball containing an NTM-Profiler database
 *
 * @output record(meta, csv, json, txt, results, logs, nf_logs, versions)
 * - `csv`: Results in CSV format
 * - `json`: Compressed machine-readable profiling results
 * - `txt`: Human-readable profiling report
 *
 * @results supplemental
 * - `*.bam`: Reads aligned to the selected species reference
 * - `*.bam.bai`: BAM alignment index
 * - `*.vcf.gz`: Variant calls and consequence annotations
 */
nextflow.enable.types = true

process NTMPROFILER_PROFILE {
    tag "${prefix}"
    label 'process_medium'

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
        csv: file("${prefix}.csv"),
        json: file("${prefix}.results.json.gz"),
        txt: file("${prefix}.txt"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.csv"),
            files("${prefix}.results.json.gz"),
            files("${prefix}.txt"),
            files("supplemental/*", optional: true)
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    def _meta = meta
    prefix = task.ext.prefix ?: "${_meta.name}"
    has_r1 = r1 != null
    has_r2 = r2 != null
    has_se = se != null
    has_lr = lr != null
    is_tarball = db.getName().endsWith(".tar.gz")
    db_setup = is_tarball ? "mkdir -p database && tar -xzf ${db} -C database" : ""
    db_arg = is_tarball ? "--db_dir database" : "--db_dir ${db}"
    input_reads = has_lr ? "--read1 ${lr}" : (has_se ? "--read1 ${se}" : "--read1 ${r1} --read2 ${r2}")
    platform = has_lr ? "--platform nanopore" : "--platform illumina"

    meta = record(
        id: "${prefix}-${task.process}",
        name: prefix,
        scope: task.ext.scope,
        output_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}",
        logs_dir: "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}",
        process_name: task.ext.process_name,
        single_end: has_se && !has_r1 && !has_r2
    )
    """
    mkdir -p results supplemental
    ${db_setup}

    ntm-profiler \\
        profile \\
        ${task.ext.args} \\
        ${platform} \\
        --csv \\
        --txt \\
        --prefix ${prefix} \\
        --threads ${task.cpus} \\
        ${db_arg} \\
        --dir results \\
        ${input_reads}

    mv results/${prefix}.results.csv ${prefix}.csv
    gzip -c results/${prefix}.results.json > ${prefix}.results.json.gz
    mv results/${prefix}.results.txt ${prefix}.txt
    for result_file in results/*.bam results/*.bai results/*.vcf.gz results/*.bcf results/*.fa; do
        if [ -e "\$result_file" ]; then
            mv "\$result_file" supplemental/
        fi
    done

    # Cleanup
    rm -rf database results

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ntm-profiler: \$(ntm-profiler --version 2>&1 | sed 's/.*version //')
    END_VERSIONS
    """
}
