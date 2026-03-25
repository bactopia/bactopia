/**
 * Profile metagenome samples against a database using Sylph.
 *
 * Uses [Sylph](https://github.com/bluenote-1/sylph) to profile metagenomic samples for taxonomic
 * abundance and containment ANI against a provided database. It is designed to be extremely fast
 * and memory-efficient.
 *
 * Uses explicit positional record fields for reads:
 * - Input: record(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords metagenomics, profiling, taxonomy, abundance, ani, sylph
 * @tags complexity:moderate input-type:single output-type:single features:conditional-logic
 * @citation sylph
 *
 * @input record(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio) - not typically used
 *
 * @input db
 * Path to the Sylph database file (*.syldb)
 *
 * @output record(meta, tsv, results, logs, nf_logs, versions)
 * - `tsv`: TSV file with profiling results
 */
nextflow.preview.types = true

process SYLPH_PROFILE {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${task.ext.container}"

    input:
    (_meta: Map, r1: Path?, r2: Path?, se: Path?, lr: Path?): Record
    db                                                      : Path

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        tsv: file("${prefix}.tsv"),
        // Generic fields (used for publishing)
        results: [
            files("${prefix}.tsv")
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
    meta.single_end = has_se && !has_r1 && !has_r2

    // Build read inputs for sylph
    def query_reads = meta.single_end ? "${se}" : "--first-pairs ${r1} --second-pairs ${r2}"
    """
    sylph \\
        profile \\
        ${db} \\
        ${query_reads} \\
        -t ${task.cpus} \\
        ${task.ext.args} \\
        --output-file ${prefix}.original.tsv

    # Remove the "fasta.gz" from sample names in output
    if [ "${meta.single_end}" == "true" ]; then
        sed 's/^${prefix}.fastq.gz/${prefix}/' ${prefix}.original.tsv > ${prefix}.tsv
    else
        sed 's/^${prefix}_R1.fastq.gz/${prefix}/' ${prefix}.original.tsv > ${prefix}.tsv
    fi
    rm ${prefix}.original.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(echo \$(sylph --version 2>&1) | sed 's/^.*sylph //;s/ .*\$//')
    END_VERSIONS
    """
}
