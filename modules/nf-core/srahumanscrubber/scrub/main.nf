// Import generic module functions
include { initOptions; saveFiles } from '../../../../lib/nf/functions'
options       = initOptions(params.containsKey("options") ? params.options : [:], 'sra-human-scrubber')
options.btype = options.btype ?: "tools"
conda_tools   = "bioconda::sra-human-scrubber=2.2.1"
conda_name    = conda_tools.replace("=", "-").replace(":", "-").replace(" ", "-")
conda_env     = file("${params.condadir}/${conda_name}").exists() ? "${params.condadir}/${conda_name}" : conda_tools

VERSION = '2.2.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
process SRAHUMANSCRUBBER_SCRUB {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-human-scrubber:2.2.1--hdfd78af_0' :
        'quay.io/biocontainers/sra-human-scrubber:2.2.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("*.scrubbed.fastq.gz"), emit: scrubbed
    path "*.{log,err}"                          , emit: logs, optional: true
    path ".command.*"                           , emit: nf_logs
    path "versions.yml"                         , emit: versions

    script:
    prefix =  options.suffix ? "${options.suffix}" : "${meta.id}"
    meta.single_end = reads[1] == null ? true : false
    if (meta.single_end) {
        """
        # Scrub human reads
        zcat ${reads} | \
            scrub.sh -d $db -p $task.cpus | \
            gzip > ${prefix}.scrubbed.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
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

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sra-human-scrubber: $VERSION
        END_VERSIONS
        """
    }
}

process SRAHUMANSCRUBBER_SCRUB_MAIN {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-human-scrubber:2.2.1--hdfd78af_0' :
        'quay.io/biocontainers/sra-human-scrubber:2.2.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads), path(extra)
    path db

    output:
    tuple val(meta), path("*.scrubbed.fastq.gz"), path(extra), emit: scrubbed
    path "*.{log,err}" , emit: logs, optional: true
    path ".command.*"  , emit: nf_logs
    path "versions.yml", emit: versions

    script:
    prefix =  options.suffix ? "${options.suffix}" : "${meta.id}"
    meta.single_end = reads[1] == null ? true : false
    if (meta.single_end) {
        """
        # Scrub human reads
        zcat ${reads} | \
            scrub.sh -d $db -p $task.cpus | \
            gzip > ${prefix}.scrubbed.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
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

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sra-human-scrubber: $VERSION
        END_VERSIONS
        """
    }
}

process SRAHUMANSCRUBBER_SCRUB_TETON {
    tag "$meta.id"
    label 'process_low'
    publishDir "${publish_dir}/${meta.id}", mode: params.publish_dir_mode, overwrite: params.force,
        saveAs: { filename -> saveFiles(filename:filename, opts:options) }
    
    conda (params.enable_conda ? conda_env : null)
    container "${ workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-human-scrubber:2.2.1--hdfd78af_0' :
        'quay.io/biocontainers/sra-human-scrubber:2.2.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads), path(extra), path(genome_size)
    path db

    output:
    tuple val(meta), path("*.scrubbed.fastq.gz"), emit: scrubbed
    path "*.{log,err}"                          , emit: logs, optional: true
    path ".command.*"                           , emit: nf_logs
    path "versions.yml"                         , emit: versions

    script:
    def prefix =  options.suffix ? "${options.suffix}" : "${meta.id}"
    meta.single_end = reads[1] == null ? true : false
    if (meta.single_end) {
        """
        # Scrub human reads
        zcat ${reads} | \
            scrub.sh -d $db -p $task.cpus | \
            gzip > ${prefix}.scrubbed.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
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

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sra-human-scrubber: $VERSION
        END_VERSIONS
        """
    }
}
