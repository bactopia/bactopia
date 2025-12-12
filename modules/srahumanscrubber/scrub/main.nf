/**
 * Scrub human reads from FASTQ files.
 *
 * Uses [SRA Human Scrubber](https://github.com/ncbi/sra-human-scrubber) to identify and remove
 * human reads from sequencing data. It relies on a specific k-mer database to mask or remove
 * sequences that align to human references.
 *
 * @status stable
 * @keywords human, contamination, scrubber, decontamination, ncbi, sra
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic
 * @citation srahumanscrubber_scrub
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: FASTQ reads to be scrubbed
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
    (_meta, reads) : Tuple<Map, List<Path>>
    db             : Path

    output:
    scrubbed             = tuple(meta, files("*.scrubbed.fastq.gz"))
    scrubbed_extra       = tuple(meta, files("*.scrubbed.fastq.gz"), file("EMPTY_EXTRA"))
    scrub_report         = tuple(meta, files('*.scrub.report.tsv', optional: true))
    scrub_special_report = tuple(special_meta, files('*.scrub.report.tsv', optional: true))
    logs                 = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs              = tuple(meta, files(".command.*"))
    versions             = tuple(meta, file("versions.yml"))

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
    meta.single_end = reads.size() == 1 ? true : false
    meta.is_paired = reads.size() == 2 ? true : false
    meta.runtype = _meta.runtype
    special_meta = [:]
    special_meta.id = prefix
    if (meta.single_end) {
        """
        # Scrub human reads
        zcat ${reads} | \
            scrub.sh -d ${db} -p ${task.cpus} | \
            gzip > ${prefix}.scrubbed.fastq.gz

        # Quick stats on reads
        zcat ${reads} | fastq-scan > original.json
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
        zcat ${reads[0]} | \
            scrub.sh -d ${db} -p ${task.cpus} | \
            gzip > ${prefix}_R1.scrubbed.fastq.gz
        zcat ${reads[1]} | \
            scrub.sh -d ${db} -p ${task.cpus} | \
            gzip > ${prefix}_R2.scrubbed.fastq.gz

        # Quick stats on reads
        zcat ${reads} | fastq-scan > original.json
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
