/**
 * k-mer based Streptococcus pneumoniae serotyping.
 *
 * Uses [SeroBA](https://github.com/sanger-pathogens/seroba) to identify the serotype of
 * *Streptococcus pneumoniae* from Illumina paired-end reads using a k-mer based approach.
 *
 * @status stable
 * @keywords streptococcus pneumoniae, serotype, k-mer, prediction, seroba
 * @tags complexity:moderate input-type:single output-type:multiple features:conditional-logic
 * @citation seroba_run
 *
 * @input tuple(meta, reads)
 * - `meta`: Groovy Map containing sample information
 * - `reads`: Paired-end FASTQ files
 *
 * @output tsv      SeroBA prediction results
 * @output txt      Detailed serogroup information
 * @output logs     Optional software execution logs containing warnings/errors
 * @output nf_logs  Nextflow execution scripts and logs for debugging
 * @output versions A YAML formatted file with software versions
 */
nextflow.preview.types = true

process SEROBA_RUN {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, reads) : Tuple<Map, Set<Path>>

    output:
    tsv      = tuple(meta, files("${prefix}.tsv"))
    txt      = tuple(meta, files("supplemental/detailed_serogroup_info.txt", optional: true))
    logs     = tuple(meta, files("*.{log,err}", optional: true))
    nf_logs  = tuple(meta, files(".command.*"))
    versions = tuple(meta, files("versions.yml"))

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
    """
    seroba \\
        runSerotyping \\
        ${reads} \\
        supplemental ${task.ext.args}

    # Avoid name collisions
    mv supplemental/pred.tsv ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seroba: \$(seroba version)
    END_VERSIONS
    """
}
