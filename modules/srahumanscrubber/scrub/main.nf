/**
 * Scrub human reads from FASTQ files.
 *
 * Uses [SRA Human Scrubber](https://github.com/ncbi/sra-human-scrubber) to identify and remove
 * human reads from sequencing data. It relies on a specific k-mer database to mask or remove
 * sequences that align to human references.
 *
 * Uses explicit positional named parameters for reads:
 * - Input: (_meta: Map, r1: Path?, r2: Path?, se: Path?, lr: Path?) as Record
 *
 * @status stable
 * @keywords human, contamination, scrubber, decontamination, ncbi, sra
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic
 * @citation srahumanscrubber
 *
 * @input (_meta: Map, r1: Path?, r2: Path?, se: Path?, lr: Path?)
 * - `_meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio) - not typically used
 *
 * @input db
 * SRA Human Scrubber database directory
 *
 * @output record(meta, special_meta, scrubbed, scrubbed_extra, scrub_report, results, logs, nf_logs, versions)
 * - `special_meta`: Groovy Map with ID for downstream aggregation
 * - `scrubbed`: Scrubbed FASTQ files with human reads removed
 * - `scrubbed_extra`: Placeholder files for pipeline compatibility
 * - `scrub_report`: Report of scrubbing statistics
 */
nextflow.preview.types = true

process SRAHUMANSCRUBBER_SCRUB {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta: Map, r1: Path?, r2: Path?, se: Path?, lr: Path?): Record
    db: Path

    output:
    record(
        // Named fields (used downstream)
        meta: meta,
        special_meta: special_meta,
        scrubbed: files("*.scrubbed.fastq.gz"),
        scrubbed_extra: files("EMPTY_EXTRA"),
        scrub_report: files('*.scrub.report.tsv', optional: true),
        // Generic fields (used for publishing)
        results: [
            files("*.scrubbed.fastq.gz"),
            files('*.scrub.report.tsv', optional: true)
        ],
        logs: files("*.{log,err}", optional: true),
        nf_logs: files(".command.*"),
        versions: files("versions.yml")
    )

    script:
    def VERSION = '2.2.1'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    prefix = task.ext.prefix ?: "${_meta.name}"
    output_folder = task.ext.wf == "scrubber" || task.ext.wf == "teton" ? "scrubber" : "${task.ext.process_name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.scope = task.ext.scope
    meta.output_dir = "${prefix}/tools/${output_folder}"
    meta.logs_dir = "${prefix}/tools/${output_folder}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name

    // Determine read type from explicit slots
    has_r1 = r1 != null
    has_r2 = r2 != null
    has_se = se != null
    meta.single_end = has_se && !has_r1 && !has_r2
    meta.runtype = _meta.runtype

    special_meta = [:]
    special_meta.id = prefix
    if (meta.single_end) {
        """
        # Scrub human reads
        zcat ${se} | \
            scrub.sh -d ${db} -p ${task.cpus} | \
            gzip > ${prefix}.scrubbed.fastq.gz

        # Quick stats on reads
        zcat ${se} | fastq-scan > original.json
        zcat *.scrubbed.fastq.gz | fastq-scan > scrubbed.json
        scrubber-summary.py ${prefix} original.json scrubbed.json > ${prefix}.scrub.report.tsv

        # Remove temp json files
        rm original.json scrubbed.json

        # Used for clean-yer-reads
        touch EMPTY_EXTRA

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastq-scan: \$(echo \$(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
            sra-human-scrubber: ${VERSION}
        END_VERSIONS
        """
    }
    else {
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
        scrubber-summary.py ${prefix} original.json scrubbed.json > ${prefix}.scrub.report.tsv

        # Remove temp json files
        rm original.json scrubbed.json

        # Used for clean-yer-reads
        touch EMPTY_EXTRA

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastq-scan: \$(echo \$(fastq-scan -v 2>&1) | sed 's/fastq-scan //')
            sra-human-scrubber: ${VERSION}
        END_VERSIONS
        """
    }
}
