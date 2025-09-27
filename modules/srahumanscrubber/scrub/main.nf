process SRAHUMANSCRUBBER_SCRUB {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.env.condaDir}/${task.ext.env.toolName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.env.image : task.ext.env.docker }"

    input:
    tuple val(_meta), path(reads)
    path db

    output:
    tuple val(meta), path("*.scrubbed.fastq.gz"), emit: scrubbed
    tuple val(meta), path("*.scrubbed.fastq.gz"), path("EMPTY_EXTRA"), emit: scrubbed_extra
    tuple val(meta), path('*.scrub.report.tsv') , emit: scrub_report, optional: true
    tuple val(meta), path("*.{log,err}")        , emit: logs, optional: true
    tuple val(meta), path(".command.begin")     , emit: nf_begin
    tuple val(meta), path(".command.err")       , emit: nf_err
    tuple val(meta), path(".command.log")       , emit: nf_log
    tuple val(meta), path(".command.out")       , emit: nf_out
    tuple val(meta), path(".command.run")       , emit: nf_run
    tuple val(meta), path(".command.sh")        , emit: nf_sh
    tuple val(meta), path(".command.trace")     , emit: nf_trace
    tuple val(meta), path("versions.yml")       , emit: versions

    script:
    def VERSION = '2.2.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    prefix = task.ext.prefix ?: "${_meta.name}"

    // Create a new meta variable
    meta = [:]
    meta.id = "${prefix}-${task.process}"
    meta.name = prefix
    meta.output_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}"
    meta.logs_dir = "${prefix}/tools/${task.ext.process_name}/${task.ext.subdir}/logs/${task.ext.logs_subdir}"
    meta.process_name = task.ext.process_name
    meta.single_end = reads[1] == null ? true : false
    meta.is_paired = reads[1] == null ? false : true
    if (meta.single_end) {
        """
        # Scrub human reads
        zcat ${reads} | \
            scrub.sh -d $db -p $task.cpus | \
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
            sra-human-scrubber: $VERSION
        END_VERSIONS
        """
    } else {
        """
        # Scrub human reads
        zcat ${reads[0]} | \
            scrub.sh -d $db -p $task.cpus | \
            gzip > ${prefix}_R1.scrubbed.fastq.gz
        zcat ${reads[1]} | \
            scrub.sh -d $db -p $task.cpus | \
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
            sra-human-scrubber: $VERSION
        END_VERSIONS
        """
    }
}
