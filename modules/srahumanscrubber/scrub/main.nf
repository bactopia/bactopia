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
    scrubbed             = tuple(meta, file("*.scrubbed.fastq.gz"))
    scrubbed_extra       = tuple(meta, file("*.scrubbed.fastq.gz"), file("EMPTY_EXTRA"))
    scrub_report         = tuple(meta, file('*.scrub.report.tsv', optional: true))
    scrub_special_report = tuple(special_meta, file('*.scrub.report.tsv', optional: true))
    logs                 = tuple(meta, file("*.{log,err}", optional: true))
    nf_begin             = tuple(meta, file(".command.begin"))
    nf_err               = tuple(meta, file(".command.err"))
    nf_log               = tuple(meta, file(".command.log"))
    nf_out               = tuple(meta, file(".command.out"))
    nf_run               = tuple(meta, file(".command.run"))
    nf_sh                = tuple(meta, file(".command.sh"))
    nf_trace             = tuple(meta, file(".command.trace"))
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
