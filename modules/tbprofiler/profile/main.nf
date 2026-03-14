/**
 * Detect resistance and lineages of Mycobacterium tuberculosis genomes.
 *
 * Uses [TBProfiler](https://github.com/jodyphelan/TBProfiler) to profile *Mycobacterium tuberculosis*
 * data for drug resistance and lineage information by aligning reads to a reference genome and identifying
 * specific variants.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords tuberculosis, mycobacterium, drug resistance, amr, typing, variant calling
 * @tags complexity:moderate input-type:single output-type:multiple features:compression,conditional-logic
 * @citation tbprofiler
 *
 * @input record(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio)
 *
 * @output record(meta, csv, json, txt, supplemental, vcf, results, logs, nf_logs, versions)
 * - `csv`: Results in CSV format
 * - `json`: Compressed JSON results file
 * - `txt`: Results in text format
 * - `supplemental`: BAM and VCF file outputs
 */
nextflow.preview.types = true

process TBPROFILER_PROFILE {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, r1: Path?, r2: Path?, se: Path?, lr: Path?): Record

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        csv: file("${prefix}.csv", optional: true),
        json: file("${prefix}.results.json.gz"),
        txt: file("${prefix}.txt", optional: true),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.csv", optional: true),
            files("${prefix}.results.json.gz"),
            files("${prefix}.txt", optional: true),
            files("supplemental/*")
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name

    // Determine read type from explicit slots
    has_r1 = r1 != null
    has_r2 = r2 != null
    has_se = se != null
    has_lr = lr != null
    meta.single_end = has_se && !has_r1 && !has_r2

    // Build read inputs and platform for tb-profiler
    def input_reads = has_lr ? "--read1 ${lr}" : (meta.single_end ? "--read1 ${se}" : "--read1 ${r1} --read2 ${r2}")
    def platform = has_lr ? "--platform nanopore" : "--platform illumina"
    """
    # Copy database to working directory
    mkdir -p database
    cp -r \$(dirname \$(which tb-profiler))/../share/tbprofiler/* database/

    tb-profiler \\
        profile \\
        ${task.ext.args} \\
        ${platform} \\
        --csv \\
        --txt \\
        --prefix ${prefix} \\
        --threads ${task.cpus} \\
        --no_trim \\
        --db_dir database/ \\
        ${input_reads}

    # Move results
    if [ -f "results/${prefix}.results.csv" ]; then
        mv results/${prefix}.results.csv ./${prefix}.csv
    fi

    # SRR2838702.results.json
    if [ -f "results/${prefix}.results.json" ]; then
        # collate hard-matches "*.results.json"
        gzip -c results/${prefix}.results.json > ${prefix}.results.json.gz
    fi

    # SRR2838702.results.txt
    if [ -f "results/${prefix}.results.txt" ]; then
        mv results/${prefix}.results.txt ./${prefix}.txt
    fi

    # Move bam and vcf folder if they exist
    mkdir supplemental
    if [ -d "bam/" ]; then
        mv bam/* supplemental/
    fi

    if [ -d "vcf/" ]; then
        mv vcf/* supplemental/
    fi

    # Cleanup
    rm -rf results/ database/ bam/ vcf/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tb-profiler:  \$( echo \$(tb-profiler profile --version 2>&1) | sed 's/.*tb-profiler version //')
    END_VERSIONS
    """
}
