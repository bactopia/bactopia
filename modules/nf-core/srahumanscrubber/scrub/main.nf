// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'sra-human-scrubber')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::bactopia-teton=1.1.0"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

VERSION = '2.2.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
process SRAHUMANSCRUBBER_SCRUB {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bactopia-teton:1.1.0--hdfd78af_0' :
        'quay.io/biocontainers/bactopia-teton:1.1.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("*.scrubbed.fastq.gz"), emit: scrubbed
    tuple val(meta), path("*.scrubbed.fastq.gz"), path("EMPTY_EXTRA"), emit: scrubbed_extra
    tuple val(meta), path('*.scrub.report.tsv') , emit: scrub_report, optional: true
    path "*.{log,err}"                          , emit: logs, optional: true
    path ".command.*"                           , emit: nf_logs
    path "versions.yml"                         , emit: versions

    script:
    prefix =  options.suffix ? "${options.suffix}" : "${meta.id}"
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
