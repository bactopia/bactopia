/**
 * Scrub human reads from FASTQ files.
 *
 * Uses [SRA Human Scrubber](https://github.com/ncbi/sra-human-scrubber) to identify and remove
 * human reads from sequencing data. It relies on a specific k-mer database to mask or remove
 * sequences that align to human references.
 *
 * Uses explicit positional tuple slots for reads:
 * - Input: tuple(meta, r1, r2, se, lr) where each read slot is Path?
 *
 * @status stable
 * @keywords human, contamination, scrubber, decontamination, ncbi, sra
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic
 * @citation srahumanscrubber
 *
 * @input tuple(meta, r1, r2, se, lr)
 * - `meta`: Groovy Map containing sample information
 * - `r1`: Illumina R1 reads (paired-end)
 * - `r2`: Illumina R2 reads (paired-end)
 * - `se`: Single-end Illumina reads
 * - `lr`: Long reads (ONT/PacBio) - not typically used
 *
 * @input db
 * SRA Human Scrubber database directory
 *
 * @output scrubbed             Scrubbed FASTQ files with human reads removed
 * @output scrubbed_extra       Scrubbed FASTQ files with placeholder for pipeline compatibility
 * @output scrub_report         Report of scrubbing statistics
 * @output scrub_special_report Special report output for downstream aggregation
 * @output logs                 Optional software execution logs containing warnings/errors
 * @output nf_logs              Nextflow execution scripts and logs for debugging
 * @output versions             A YAML formatted file with software versions
 */
nextflow.preview.types = true

process SRAHUMANSCRUBBER_SCRUB {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, r1, r2, se, lr) : Tuple<Map, Path?, Path?, Path?, Path?>
    db                      : Path

    output:
    scrubbed             = tuple(meta, files("*.scrubbed.fastq.gz"))
    scrubbed_extra       = tuple(meta, files("*.scrubbed.fastq.gz"), files("EMPTY_EXTRA"))
    scrub_report         = tuple(meta, files('*.scrub.report.tsv', optional: true))
    scrub_special_report = tuple(special_meta, files('*.scrub.report.tsv', optional: true))
    logs                 = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs              = tuple(meta, files(".command.*"))
    versions             = tuple(meta, files("versions.yml"))

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
