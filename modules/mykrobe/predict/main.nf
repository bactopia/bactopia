/**
 * Predict Antimicrobial Resistance (AMR) for supported bacterial species.
 *
 * Uses [Mykrobe](https://github.com/mykrobe/mykrobe) to quickly predict resistance and susceptibility
 * based on short reads (FASTQ) or aligned reads (BAM). It maps k-mers from the input sequences
 * against a curated database of resistance markers for species like *M. tuberculosis* and *S. aureus*.
 *
 * Uses explicit positional record fields for reads:
 * - Input: record(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords amr, resistance, susceptibility, k-mer, fastq, bam, mykrobe
 * @tags complexity:moderate input-type:multiple output-type:multiple features:database-dependent
 * @citation mykrobe, mccortex
 *
 * @input record(meta, r1?, r2?, se?, lr?)
 * - `meta`: Groovy Record containing sample information
 * - `r1?`: Illumina R1 reads (paired-end)
 * - `r2?`: Illumina R2 reads (paired-end)
 * - `se?`: Single-end Illumina reads
 * - `lr?`: Long reads (ONT/PacBio)
 *
 * @input species
 * The target species for which to make the AMR prediction (e.g., "tb" or "staph")
 *
 * @output record(meta, csv, json, results, logs, nf_logs, versions)
 * - `csv`: AMR predictions in machine-readable CSV format
 * - `json`: Detailed AMR prediction results in JSON format
 */
nextflow.enable.types = true

process MYKROBE_PREDICT {
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
    species: String

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        csv: file("${prefix}.csv"),
        json: file("${prefix}.json"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.csv"),
            files("${prefix}.json")
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
        single_end: has_se && !has_r1 && !has_r2
    )

    // Build read inputs and options for mykrobe
    def is_ont = has_lr ? "--ont" : ""
    def read_inputs = has_lr ? "${lr}" : (meta.single_end ? "${se}" : "${r1} ${r2}")
    """
    mykrobe \\
        predict \\
        ${task.ext.args} ${is_ont} \\
        --species ${species} \\
        --threads ${task.cpus} \\
        --sample ${prefix} \\
        --format json_and_csv \\
        --output ${prefix} \\
        --seq ${read_inputs}

    # Cleanup
    rm -rf mykrobe

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mykrobe: \$(echo \$(mykrobe --version 2>&1) | sed 's/^.*mykrobe v//' )
    END_VERSIONS
    """
}
