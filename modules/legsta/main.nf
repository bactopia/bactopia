/**
 * In silico Sequence Based Typing (SBT) of *Legionella pneumophila*.
 *
 * Uses [Legsta](https://github.com/tseemann/legsta) to determine the Sequence Based Type (SBT)
 * of *L. pneumophila* isolates. It aligns the assembly against the standard 7-gene schema
 * (flaA, pilE, asd, mip, mompS, proA, neuA) to assign allele numbers and the resulting Sequence Type.
 *
 * @status stable
 * @keywords bacteria, legionella, pneumophila, typing, sbt, mlst, serogroup
 * @tags complexity:simple input-type:single output-type:single
 * @citation legsta
 *
 * @input tuple(meta, assembly)
 * - `meta`: Groovy Map containing sample information
 * - `assembly`: Assembled contigs in FASTA format
 *
 * @output tsv       A tab-delimited summary of the assigned Sequence Type and allele profiles
 * @output logs      Optional software execution logs containing warnings/errors
 * @output nf_logs   Nextflow execution scripts and logs for debugging
 * @output versions  A YAML formatted file with software versions
 */
nextflow.preview.types = true

process LEGSTA {
    tag "${prefix}"
    label 'process_low'

    conda "${task.ext.condaDir}/${task.ext.toolName}"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? task.ext.image : task.ext.docker}"

    input:
    (_meta, assembly) : Tuple<Map, Path>

    output:
    tsv      = tuple(meta, file("${prefix}.tsv"))
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
    legsta \\
        ${task.ext.args} \\
        ${assembly} | sed 's/.fna//; s/.gz//' > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        legsta: \$(echo \$(legsta --version 2>&1) | sed 's/^.*legsta //; s/ .*\$//;')
    END_VERSIONS
    """
}
