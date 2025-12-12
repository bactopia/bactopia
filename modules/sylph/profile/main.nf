/**
 * Profile metagenome samples against a database using Sylph.
 *
 * Uses [Sylph](https://github.com/bluenote-1/sylph) to profile metagenomic samples for taxonomic
 * abundance and containment ANI against a provided database. It is designed to be extremely fast
 * and memory-efficient.
 *
 * @status stable
 * @keywords metagenomics, profiling, taxonomy, abundance, ani, sylph
 * @tags complexity:moderate input-type:single output-type:single features:conditional-logic
 * @citation sylph_profile
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: FASTQ reads (single or paired-end)
 *
 * @input db
 * Path to the Sylph database file (*.syldb)
 *
 * @output tsv      TSV file with profiling results
 * @output logs     Optional software execution logs containing warnings/errors
 * @output nf_logs  Nextflow execution scripts and logs for debugging
 * @output versions A YAML formatted file with software versions
 */
nextflow.preview.types = true

process SYLPH_PROFILE {
    tag "${prefix}"
    label 'process_medium'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, reads) : Tuple<Map, List<Path>>
    db             : Path

    output:
    tsv      = tuple(meta, file("${prefix}.tsv"))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, file("versions.yml"))

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
    def query_reads = meta.single_end ? "${reads[0]}" : "--first-pairs ${reads[0]} --second-pairs ${reads[1]}"
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
